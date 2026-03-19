/******************************************************************************
 *   Copyright (C) 2014 - 2020 Jan Fostier (jan.fostier@ugent.be)             *
 *   This file is part of Detox                                               *
 *                                                                            *
 *   This program is free software; you can redistribute it and/or modify     *
 *   it under the terms of the GNU Affero General Public License as published *
 *   by the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                      *
 *                                                                            *
 *   This program is distributed in the hope that it will be useful,          *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *   GNU General Public License for more details.                             *
 *                                                                            *
 *   You should have received a copy of the GNU Affero General Public License *
 *   along with this program; if not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************/

#include <numeric>
#include <thread>

#include "crfmult.h"
#include "dbgraph.h"
#include "pgm/pgminference.h"

using namespace std;

double CRFSolver::logFS = 16.118095651;  // log(1e7)

// ============================================================================
// CONDITIONAL RANDOM FIELDS SOLVER (PER THREAD)
// ============================================================================

void CRFSolver::getSubgraph(int maxDepth)
{
        while (!pq.empty()) {
                // get and erase the current node
                auto [thisID, thisDepth] = pq.top();
                pq.pop();

                // if the node was already handled, skip
                if (bitvec[abs(thisID)])
                        continue;

                nodes.push_back(NodeRep(thisID));
                SSNode n = dBG.getSSNode(thisID);

                // process the right arcs
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        NodeID rightID = it->getNodeID();
                        if (bitvec[abs(rightID)])       // edge already added?
                                continue;

                        edges.push_back(EdgeRep(thisID, rightID));
                        if (thisDepth < maxDepth)
                                pq.push(NodeDepth(rightID, thisDepth+1));
                }

                // process the left arcs
                for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++) {
                        NodeID leftID = it->getNodeID();
                        if (bitvec[abs(leftID)])        // edge already added?
                                continue;
                        if (leftID == thisID)   // edge already added as right arc
                                continue;

                        edges.push_back(EdgeRep(leftID, thisID));
                        if (thisDepth < maxDepth)
                                pq.push(NodeDepth(leftID, thisDepth + 1));
                }

                // mark this node as handled
                bitvec[abs(thisID)] = true;
        }

        // reset all flags to false
        for (auto it : nodes)
                bitvec[it.getNodeID()] = false;
}

void CRFSolver::getSubgraph(NodeRep seedNode, int maxDepth)
{
        nodes.clear(); edges.clear();

        if (maxDepth == 0) {
                nodes.push_back(seedNode);
                return;
        }

        // add the seed node to the priority queue
        pq.push(NodeDepth(seedNode, 0));

        getSubgraph(maxDepth);
}

void CRFSolver::getSubgraph(EdgeRep seedEdge, int maxDepth)
{
        nodes.clear(); edges.clear();

        // special case (maxDepth == 0)
        if (maxDepth == 0) {
                edges.push_back(seedEdge);
                return;
        }

        // other cases (maxDepth > 0)
        pq.push(NodeDepth(seedEdge.getSrcID(), 1));
        pq.push(NodeDepth(seedEdge.getDstID(), 1));

        getSubgraph(maxDepth);
}

Factor CRFSolver::createSingletonFactor(int varID, int varCard, int firstMult,
                                        const CovModel& covModel, double obsCov,
                                        int numIndepObs, bool tip)
{
        vector<int> var(1, varID);
        vector<int> card(1, varCard);
        vector<double> val(varCard);

        // create a singleton
        for (int i = 0; i < varCard; i++) {
                int multiplicity = firstMult + i;
                val[i] = covModel.getLogProb(obsCov, multiplicity);
        }

        // F = F^numIndepObs (in log-space)
        for (int i = 0; i < varCard; i++)
                val[i] *= numIndepObs;

        // express that flow is not conserved at tips
        if (tip && (firstMult == 0))
                for (int i = 1; i < varCard; i++)
                        val[i] += getLogFlowStrength(0, i);

        return Factor(move(var), move(card), move(val));
}

Factor CRFSolver::createFlowFactor(int sumVarID, int sumCard, int sumFirstMult,
                                   const vector<int>& termVarID,
                                   const vector<int>& termCard,
                                   const vector<int>& termFirstMult,
                                   const vector<bool>& palindromic)
{
        assert(termVarID.size() == termFirstMult.size());

        // create the variable vector
        vector<int> var = { sumVarID };
        var.insert(var.end(), termVarID.begin(), termVarID.end());

        // create the cardinality vector
        vector<int> card = { sumCard };
        card.insert(card.end(), termCard.begin(), termCard.end());

        // create the value vector
        size_t numVal = accumulate(card.begin(), card.end(), 1ull, multiplies<size_t>());
        vector<double> val(numVal);

        Assignment assignment(card);
        for (size_t i = 0; i < val.size(); i++) {
                int centralVal = assignment[0] + sumFirstMult;
                int sumTerms = 0;
                for (size_t j = 0; j < termFirstMult.size(); j++) {
                        // count palindromic arcs double
                        int c = (palindromic[j]) ? 2 : 1;
                        sumTerms += c * (termFirstMult[j] + assignment[j+1]);
                }

                val[i] = getLogFlowStrength(centralVal, sumTerms);
                assignment.increment();
        }

        return Factor(move(var), move(card), move(val));
}

bool CRFSolver::solveSubgraph(const map<NodeRep, int>& node2var,
                              const map<EdgeRep, int>& edge2var,
                              const CovModel& nodeCovModel,
                              const CovModel& edgeCovModel,
                              Multiplicity& resMult, int targetCRFVar) const
{
        // create the singleton factors
        list<Factor> factorList;
        vector<int> card(node2var.size() + edge2var.size());
        vector<int> firstMult(node2var.size() + edge2var.size(), -1);

        for (const auto [nr, nodeVarID] : node2var) {   // first the nodes
                // always compute firstMult and card even if factor is not added
                SSNode node = dBG.getSSNode(nr);
                int expMult = nodeCovModel.getExpMult(node.getAvgCov());
                firstMult[nodeVarID] = max(0, expMult - multMargin);
                card[nodeVarID] = expMult + multMargin - firstMult[nodeVarID] + 1;

                //cout << "Node: " << nr << ", CRF var: " << nodeVarID << endl;

                // add a singleton node factor in case:
                // - the node is at least 2k (less correlation with edges) OR
                // - there are no edges (in case of isolated nodes) OR
                // - there is a single incoming edge OR
                // - there is a single outgoing edge
                if ((node.getMarginalLength() < 2*Kmer::getK()) &&
                    (node.numLeftArcs() > 1) && (node.numRightArcs() > 1))
                        continue;

                // the number of indep observations accounts for node length
                const NodeLength& nl = node.getMarginalLength();
                int nIndptObs = max<int>(1, nl / INDEP_OBS);

                // at a tip, flow is not conserved. This is accounted for in
                // the singleton factor itself.
                bool tip = ((node.numLeftArcs() == 0) ||
                            (node.numRightArcs() == 0));
                tip = tip && (node.getMarginalLength() <= MAX_TIP_LEN);

                //cout << "Singleton NODE added" << endl;
                Factor F = createSingletonFactor(nodeVarID, card[nodeVarID],
                                                 firstMult[nodeVarID],
                                                 nodeCovModel, node.getAvgCov(),
                                                 nIndptObs, tip);
                factorList.push_back(F);
        }

        for (const auto [er, edgeVarID] : edge2var) {   // then the edges
                // compute firstMult and card even if factor is not added
                double cov = dBG.getArc(er).getCov();

                if (firstMult[edgeVarID] == -1) {
                        int expMult = edgeCovModel.getExpMult(cov);
                        firstMult[edgeVarID] = max(0, expMult - multMargin);
                        card[edgeVarID] = expMult + multMargin - firstMult[edgeVarID] + 1;
                }

                //cout << "Edge: " << er << ", var: " << edgeVarID << endl;

                // add a singleton edge factor in case:
                // - an adjacent node is not in the CRF (extremal edge) OR
                // - both adjacent nodes have multiple edges on that side
                if (node2var.find(NodeRep(er.getSrcID())) != node2var.end()) {
                        SSNode src = dBG.getSSNode(er.getSrcID());
                        if (src.numRightArcs() == 1)
                                continue;
                }
                if (node2var.find(NodeRep(er.getDstID())) != node2var.end()) {
                        SSNode dst = dBG.getSSNode(er.getDstID());
                        if (dst.numLeftArcs() == 1)
                                continue;
                }

                //cout << "Singleton EDGE added" << endl;
                Factor F = createSingletonFactor(edgeVarID, card[edgeVarID],
                                                 firstMult[edgeVarID],
                                                 edgeCovModel, cov, 1, false);
                factorList.push_back(F);
        }

        // create the flow conservation factors
        if ( !edge2var.empty() ){
        for (const auto it : node2var) {
                NodeID currID = it.first.getNodeID();
                const int& nodeVarID = it.second;
                SSNode node = dBG.getSSNode(currID);

                // left flow conservation factor
                vector<int> lEdgeVarID, lEdgeVarCard, lEdgeFirstMult;
                vector<bool> lEdgePalindromic;
                for (ArcIt lIt = node.leftBegin(); lIt != node.leftEnd(); lIt++) {
                        EdgeRep edge(lIt->getNodeID(), currID);
                        int edgeVarID = edge2var.at(edge);

                        // Node has a single edge at that side with the same
                        // CRF variable. No need for flow conservation factor.
                        if (edgeVarID == nodeVarID)
                                continue;

                        lEdgeVarID.push_back(edgeVarID);
                        lEdgeVarCard.push_back(card[edgeVarID]);
                        lEdgeFirstMult.push_back(firstMult[edgeVarID]);
                        lEdgePalindromic.push_back(lIt->getNodeID() == -currID);
                }

                if (!lEdgeVarID.empty()) {
                        Factor F = createFlowFactor(nodeVarID, card[nodeVarID],
                                                    firstMult[nodeVarID], lEdgeVarID,
                                                    lEdgeVarCard, lEdgeFirstMult,
                                                    lEdgePalindromic);
                        factorList.push_back(F);
                }

                // right flow conservation factor
                vector<int> rEdgeVarID, rEdgeVarCard, rEdgeFirstMult;
                vector<bool> rEdgePalindromic;
                for (ArcIt rIt = node.rightBegin(); rIt != node.rightEnd(); rIt++) {
                        EdgeRep edge(currID, rIt->getNodeID());
                        int edgeVarID = edge2var.at(edge);

                        // Node has a single edge at that side with the same
                        // CRF variable. No need for flow conservation factor.
                        if (edgeVarID == nodeVarID)
                                continue;

                        rEdgeVarID.push_back(edgeVarID);
                        rEdgeVarCard.push_back(card[edgeVarID]);
                        rEdgeFirstMult.push_back(firstMult[edgeVarID]);
                        rEdgePalindromic.push_back(currID == -rIt->getNodeID());
                }

                if (!rEdgeVarID.empty()) {
                        Factor F = createFlowFactor(nodeVarID, card[nodeVarID],
                                                    firstMult[nodeVarID], rEdgeVarID,
                                                    rEdgeVarCard, rEdgeFirstMult,
                                                    rEdgePalindromic);
                        factorList.push_back(F);
                }
        }
        }

        /*for (auto it : factorList) {
                cout << it << endl;
                cout << " ==================== " << endl;
        }*/

        // Variable elimination
        bool retVal = PGMInference::solveVE(factorList, targetCRFVar, maxFactorSize);
        if (!retVal)
                return false;

        resMult = Multiplicity(firstMult[targetCRFVar],
                               factorList.front().getVal());
        resMult.normalize();

        //cout << "RESULT: " << resMult << endl;

        return true;
}

bool CRFSolver::checkFlow(NodeRep node, const CovModel& nodeCovModel,
                          const CovModel& edgeCovModel, int graphDepth)
{
        if (!dBG.getSSNode(node).isValid())
                return false;
        getSubgraph(node, graphDepth);

        for (auto e : nodes) {
                SSNode n = dBG.getSSNode(e);
                int nodeMult = nodeCovModel.getExpMult(n.getAvgCov());

                int sumRight = 0;
                for (ArcIt r = n.rightBegin(); r != n.rightEnd(); r++) {
                        // count palindromic arcs double
                        int c = (r->getNodeID() == -e) ? 2 : 1;
                        sumRight += c * edgeCovModel.getExpMult(r->getCov());
                }

                if ((n.numRightArcs() > 0) && (nodeMult != sumRight))
                        return false;

                int sumLeft = 0;
                for (ArcIt l = n.leftBegin(); l != n.leftEnd(); l++) {
                        // count palindromic arcs double
                        int c = (l->getNodeID() == -e) ? 2 : 1;
                        sumLeft += c * edgeCovModel.getExpMult(l->getCov());
                }

                if ((n.numLeftArcs() > 0) && (nodeMult != sumLeft))
                        return false;
        }

        return true;
}

bool CRFSolver::checkFlow(EdgeRep edge, const CovModel& nodeCovModel,
                          const CovModel& edgeCovModel, int graphDepth)
{
        getSubgraph(edge, graphDepth);

        for (auto e : nodes) {
                SSNode n = dBG.getSSNode(e);
                int nodeMult = nodeCovModel.getExpMult(n.getAvgCov());

                int sumRight = 0;
                for (ArcIt r = n.rightBegin(); r != n.rightEnd(); r++) {
                        // count palindromic arcs double
                        int c = (r->getNodeID() == -e) ? 2 : 1;
                        sumRight += c * edgeCovModel.getExpMult(r->getCov());
                }

                if ((n.numRightArcs() > 0) && (nodeMult != sumRight))
                        return false;

                int sumLeft = 0;
                for (ArcIt l = n.leftBegin(); l != n.leftEnd(); l++) {
                        // count palindromic arcs double
                        int c = (l->getNodeID() == -e) ? 2 : 1;
                        sumLeft += c * edgeCovModel.getExpMult(l->getCov());
                }

                if ((n.numLeftArcs() > 0) && (nodeMult != sumLeft))
                        return false;
        }

        return true;
}

Multiplicity CRFSolver::getNodeMultiplicity(NodeRep node,
                                            const CovModel& nodeCovModel,
                                            const CovModel& edgeCovModel,
                                            int graphDepth)
{
        getSubgraph(node, graphDepth);

        // create a mapping between nodes/edges and CRF variables
        map<NodeRep, int> node2var;
        map<EdgeRep, int> edge2var;

        // Each node is assigned a unique varCRF. Edges incident to nodes
        // that have (on that side) only single edge are assigned the same
        // CRF variable as that node. This means that a node, together with
        // its unique edge (or edges: in case of a unique incoming and
        // outgoing edge) will always have the same multiplicity.
        int varCRF = 0;
        for (const auto it : nodes) {
                node2var[it] = varCRF;
                SSNode n = dBG.getSSNode(it);
                if (n.numLeftArcs() == 1) {     // edge gets same varCRF
                        EdgeRep er(n.leftBegin()->getNodeID(), it);
                        edge2var[er] = varCRF;
                }
                if (n.numRightArcs() == 1) {    // edge gets same varCRF
                        EdgeRep er(it, n.rightBegin()->getNodeID());
                        edge2var[er] = varCRF;
                }
                varCRF++;
        }

        // Assign a varCRF to edges that are still unassigned
        for (const auto it : edges) {
                if (edge2var.find(it) == edge2var.end())
                        edge2var[it] = varCRF++;
        }

        Multiplicity result;
        int targetCRFVar = node2var[node];
        if (solveSubgraph(node2var, edge2var, nodeCovModel,
                          edgeCovModel, result, targetCRFVar)) {
                perfCounter.totNumCRF++;
                perfCounter.totNumNodes += nodes.size();
                perfCounter.totNumEdges += edges.size();

                return result;
        }

        // fall-back to smaller subgraph if necessary
        return getNodeMultiplicity(node, nodeCovModel,
                                   edgeCovModel, graphDepth - 1);
}

Multiplicity CRFSolver::getEdgeMultiplicity(EdgeRep edge,
                                            const CovModel& nodeCovModel,
                                            const CovModel& edgeCovModel,
                                            int graphDepth)
{
        getSubgraph(edge, graphDepth);

        // create a mapping between nodes/edges and CRF variables
        map<NodeRep, int> node2var;
        map<EdgeRep, int> edge2var;

        // Each node is assigned a unique varCRF. Edges incident to nodes
        // that have (on that side) only single edge are assigned the same
        // CRF variable as that node. This means that a node, together with
        // its unique edge (or edges: in case of a unique incoming and
        // outgoing edge) will always have the same multiplicity.
        int varCRF = 0;
        for (const auto it : nodes) {
                node2var[it] = varCRF;
                SSNode n = dBG.getSSNode(it);
                if (n.numLeftArcs() == 1) {     // edge gets same varCRF
                        EdgeRep er(n.leftBegin()->getNodeID(), it);
                        edge2var[er] = varCRF;
                }
                if (n.numRightArcs() == 1) {    // edge gets same varCRF
                        EdgeRep er(it, n.rightBegin()->getNodeID());
                        edge2var[er] = varCRF;
                }
                varCRF++;
        }

        // Assign a varCRF to edges that are still unassigned
        for (const auto it : edges) {
                if (edge2var.find(it) == edge2var.end())
                        edge2var[it] = varCRF++;
        }

        Multiplicity result;
        int targetCRFVar = edge2var[edge];
        if (solveSubgraph(node2var, edge2var, nodeCovModel,
                          edgeCovModel, result, targetCRFVar))
                return result;

        // fall-back to smaller subgraph if necessary
        return getEdgeMultiplicity(edge, nodeCovModel, edgeCovModel, graphDepth-1);
}

// ============================================================================
// ESTIMATE NODE/EDGE MULTIPLICITIES USING CONDITIONAL RANDOM FIELDS
// ============================================================================

void CRFMult::checkNodeFlow(CRFSolver& solver, WorkLoadBalancer& wlb,
                            const vector<NodeRep>& nodes,
                            const CovModel& nodeCovModel,
                            const CovModel& edgeCovModel,
                            vector<bool>& flowOK) const
{
        size_t idxBegin, idxEnd;
        while (wlb.getChunk(idxBegin, idxEnd))
                for (size_t id = idxBegin; id < idxEnd; id++)
                        flowOK[id] = solver.checkFlow(nodes[id],
                                                      nodeCovModel,
                                                      edgeCovModel,
                                                      maxGraphDepth);
}

void CRFMult::checkEdgeFlow(CRFSolver& solver, WorkLoadBalancer& wlb,
                           const std::vector<EdgeRep>& edges,
                           const CovModel& nodeCovModel,
                           const CovModel& edgeCovModel,
                           std::vector<bool>& flowOK) const
{
        size_t idxBegin, idxEnd;
        while (wlb.getChunk(idxBegin, idxEnd))
                for (size_t id = idxBegin; id < idxEnd; id++)
                        flowOK[id] = solver.checkFlow(edges[id],
                                                      nodeCovModel,
                                                      edgeCovModel,
                                                      maxGraphDepth);
}

void CRFMult::computeNodeMult(CRFSolver& solver,
                              WorkLoadBalancer& wlb,
                              const vector<NodeRep>& nodes,
                              const CovModel& nodeCovModel,
                              const CovModel& edgeCovModel,
                              NodeMap<Multiplicity>& nodeMult) const
{
        size_t idxBegin, idxEnd;
        while (wlb.getChunk(idxBegin, idxEnd))
                for (size_t id = idxBegin; id < idxEnd; id++) {
                        NodeRep n = nodes[id];
                        nodeMult[n] = solver.getNodeMultiplicity(n,
                                nodeCovModel, edgeCovModel, maxGraphDepth);
                }
}

void CRFMult::computeEdgeMult(CRFSolver& solver,
                              WorkLoadBalancer& wlb,
                              const vector<EdgeRep>& edges,
                              const CovModel& nodeCovModel,
                              const CovModel& edgeCovModel,
                              EdgeMap<Multiplicity>& edgeMult) const
{
        size_t idxBegin, idxEnd;
        while (wlb.getChunk(idxBegin, idxEnd))
                for (size_t id = idxBegin; id < idxEnd; id++) {
                        const EdgeRep& e = edges[id];
                        edgeMult[e] = solver.getEdgeMultiplicity(e,
                                nodeCovModel, edgeCovModel, maxGraphDepth);
                }
}

void CRFMult::checkFlow(const vector<NodeRep>& nodes,
                        vector<bool>& flowOK,
                        const CovModel& nodeCovModel,
                        const CovModel& edgeCovModel) const
{
        cout << "Checking flow for " << nodes.size() << " nodes (subgraph "
                "depth: " << maxGraphDepth << ")" << endl;
        flowOK.resize(nodes.size());
        vector<thread> wt(numThreads);
        vector<CRFSolver> solver(numThreads, CRFSolver(dBG, multMargin,
                                                       maxFactorSize, flowStrength));

        // assign a multiplicity to the nodes
        WorkLoadBalancer nodeWlb(0, nodes.size(), threadWork, "\tProcessing nodes");
        for (size_t i = 0; i < wt.size(); i++)
                wt[i] = thread(&CRFMult::checkNodeFlow, this, ref(solver[i]),
                               ref(nodeWlb), cref(nodes), cref(nodeCovModel),
                               cref(edgeCovModel), ref(flowOK));

        // wait for worker threads to finish
        for_each(wt.begin(), wt.end(), mem_fn(&thread::join));
}

void CRFMult::checkFlow(const vector<EdgeRep>& edges,
                        vector<bool>& flowOK,
                        const CovModel& nodeCovModel,
                        const CovModel& edgeCovModel) const
{
        cout << "Checking flow for " << edges.size() << " edges (subgraph "
                "depth: " << maxGraphDepth << ")" << endl;
        flowOK.resize(edges.size());
        vector<thread> wt(numThreads);
        vector<CRFSolver> solver(numThreads, CRFSolver(dBG, multMargin,
                                                       maxFactorSize, flowStrength));

        // assign a multiplicity to the nodes
        WorkLoadBalancer nodeWlb(0, edges.size(), threadWork, "\tProcessing edges");
        for (size_t i = 0; i < wt.size(); i++)
                wt[i] = thread(&CRFMult::checkEdgeFlow, this, ref(solver[i]),
                               ref(nodeWlb), cref(edges), cref(nodeCovModel),
                               cref(edgeCovModel), ref(flowOK));

        // wait for worker threads to finish
        for_each(wt.begin(), wt.end(), mem_fn(&thread::join));
}

void CRFMult::computeMult(NodeMap<Multiplicity>& nodeMult,
                          const CovModel& nodeCovModel,
                          EdgeMap<Multiplicity>& edgeMult,
                          const CovModel& edgeCovModel, bool silent) const
{
        if (!silent)
                cout << "Computing multiplicity for " << nodeMult.size()
                     << " nodes and " << edgeMult.size() << " edges (subgraph "
                        "depth: " << maxGraphDepth << ")" << endl;
        vector<thread> wt(numThreads);
        vector<CRFSolver> solver(numThreads,
                CRFSolver(dBG, multMargin, maxFactorSize, flowStrength));

        // assign a multiplicity to the nodes
        vector<NodeRep> nodes; nodes.reserve(nodeMult.size());
        for (const auto& it : nodeMult)
                nodes.push_back(it.first);
        string message = silent ? string() : "\tProcessing nodes";
        WorkLoadBalancer nodeWlb(0, nodes.size(), threadWork, message);
        for (size_t i = 0; i < wt.size(); i++)
                wt[i] = thread(&CRFMult::computeNodeMult, this, ref(solver[i]),
                               ref(nodeWlb), cref(nodes), cref(nodeCovModel),
                               cref(edgeCovModel), ref(nodeMult));

        // wait for worker threads to finish
        for_each(wt.begin(), wt.end(), mem_fn(&thread::join));
        nodes.clear();

        // assign a multiplicity to the edges
        vector<EdgeRep> edges; edges.reserve(edgeMult.size());
        for (const auto& it : edgeMult)
                edges.push_back(it.first);
        message = silent ? string() : "\tProcessing edges";
        WorkLoadBalancer edgeWlb(0, edges.size(), threadWork, message);
        for (size_t i = 0; i < wt.size(); i++)
                wt[i] = thread(&CRFMult::computeEdgeMult, this, ref(solver[i]),
                               ref(edgeWlb), cref(edges), cref(nodeCovModel),
                               cref(edgeCovModel), ref(edgeMult));

        // wait for worker threads to finish
        for_each(wt.begin(), wt.end(), mem_fn(&thread::join));

        // merge the performance counters
        for (const auto& it : solver)
                totPerfCounter.merge(it.perfCounter);
}

bool CRFMult::Estep(NodeMap<Multiplicity>& nodeMult,
                    const CovModel& nodeCovModel,
                    EdgeMap<Multiplicity>& edgeMult,
                    const CovModel& edgeCovModel, double epsilon) const
{
        // store the old expected multiplicities to check for convergence
        NodeMap<int> oldNodeMult(nodeMult.size());
        for (const auto& it : nodeMult)
                oldNodeMult[it.first] = it.second.getExpMult();
        EdgeMap<int> oldEdgeMult(edgeMult.size());
        for (const auto& it : edgeMult)
                oldEdgeMult[it.first] = it.second.getExpMult();

        // compute (normalized) multiplicities given the model
        computeMult(nodeMult, nodeCovModel, edgeMult, edgeCovModel);

        // check for convergence
        size_t numNodeChange = 0;
        for (auto it : nodeMult)
                if (oldNodeMult[it.first] <= 5)
                if (it.second.getExpMult() != oldNodeMult[it.first])
                        numNodeChange++;
        size_t numEdgeChange = 0;
        for (auto it : edgeMult)
                if (oldEdgeMult[it.first] <= 5)
                if (it.second.getExpMult() != oldEdgeMult[it.first])
                        numEdgeChange++;

        // convergence
        return (((double)numNodeChange/nodeMult.size() > epsilon) ||
                ((double)numEdgeChange/edgeMult.size() > epsilon));
}

bool CRFMult::Mstep(const NodeMap<Multiplicity>& nodeMult,
                    CovModel& nodeCovModel,
                    const EdgeMap<Multiplicity>& edgeMult,
                    CovModel& edgeCovModel, bool trainNodeODF,
                    double epsilon) const
{
        // ====================================================================
        // Compute the node model parameters
        // ====================================================================

        // We fit a negative binomial (NB) distribution to the error histogram.
        // The error histogram is always truncated: a) coverage zero is never
        // observed and b) low-coverage k-mers might have been removed by
        // a preprocessing tool like BCALM. We therefore use EM to fit the NB
        // parameters and infer the missing values from the spectrum.
        map<unsigned int, double> errorHist;
        for (const auto& it : nodeMult) {
                SSNode node = dBG.getSSNode(it.first);
                if (it.second[0] > DOUBLE_SMALL) {
                        double f = node.getAvgCov() - floor(node.getAvgCov());
                        errorHist[node.getAvgCov()] += (1.0-f) * it.second[0];
                        errorHist[node.getAvgCov() + 1] += f * it.second[0];
                }
        }

        // We truncate the node histogram to the provided -abundance-min value.
        // This is relevant when using q-mer counts.
        unsigned int smallSupport = (dBG.getAbundanceMin() < 0)?
                errorHist.begin()->first : dBG.getAbundanceMin();
        for (unsigned int k = 0; k < smallSupport; k++)
                errorHist.erase(k);

        double nodeErrorLambda = nodeCovModel.getErrLambda();
        double nodeErrorODF = nodeCovModel.getErrorODF();
        double nodeErrorWeight = nodeCovModel.getWeight(0);

        // The below EM procedure might not convergence, but the obtained NB
        // parameters should be sensible in all cases: only small support
        // values are inferred, therefore, the mean should never be larger than
        // the mean of the error histogram (and thus be close to zero).
        int maxIter = 10.0 / epsilon;
        int nIter = Util::fitTruncNegBinomEM(errorHist, nodeErrorLambda,
                                             nodeErrorODF, nodeErrorWeight,
                                             epsilon, maxIter);

        if (nIter > maxIter)
                cout << "\tWARNING: EM algorithm to fit node error model "
                        "did not converge\n";
        else
                cout << "\tEM algorithm to fit node error model converged "
                     << "after " << nIter << " iteration(s)" << endl;

        // Compute node average and weights
        vector<double> wNode(nodeCovModel.getK(), 0.0);
        wNode[0] = nodeErrorWeight;

        double noml = 0.0, denl = 0.0;
        for (const auto& it : nodeMult) {
                SSNode node = dBG.getSSNode(it.first);
                for (size_t m = 1; m < nodeCovModel.getK(); m++) {
                        wNode[m] += it.second[m];
                        noml += it.second[m] * node.getCov();
                        denl += it.second[m] * node.getMarginalLength() * m;
                }
        }

        double nodeLambda = max(0.5, noml / denl);

        // Compute node variance and ODF
        vector<double> wNodeSqr(nodeCovModel.getK(), 0.0);

        for (const auto& it : nodeMult) {
                SSNode node = dBG.getSSNode(it.first);
                for (size_t m = 1; m < nodeCovModel.getK(); m++) {
                        double delta = node.getAvgCov() - m * nodeLambda;
                        wNodeSqr[m] += it.second[m] * delta * delta;
                }
        }
        double nom = 0.0, den = 0.0;
        for (size_t m = 1; m < nodeCovModel.getK(); m++) {
                nom += wNodeSqr[m] / (m * nodeLambda);
                den += wNode[m];
        }
        double nodeODF = nom / den;
        nodeODF = max(1.0, nodeODF);

        // The node ODF should not be retrained if the graph contains many
        // long nodes. The avgCov() of those nodes will already be averaged
        // over many k-mers and therefore, the ODF will be underestimated.
        // In turn, this will lead to too confident singleton factors.
        if (!trainNodeODF)
                nodeODF = nodeCovModel.getODF();

        nodeCovModel = CovModel(nodeErrorLambda, nodeErrorODF,
                                nodeLambda, nodeODF, wNode);

        cout << "\tNode spectrum: " << nodeCovModel << endl;

        // ====================================================================
        // Compute the edge model parameters
        // ====================================================================

        errorHist.clear();
        for (const auto& it : edgeMult) {
                Arc& arc = dBG.getArc(it.first);
                if (it.second[0] > DOUBLE_SMALL) {
                        double f = arc.getCov() - floor(arc.getCov());
                        errorHist[arc.getCov()] += (1.0-f) * it.second[0];
                        errorHist[arc.getCov() + 1] += f * it.second[0];
                }
        }

        // We truncate the edge histogram on the same value as the node
        // histogram. Even when removing all k-mers with coverage < T, a few
        // arcs might still have a coverage < T. We do not want keep those
        // in our histogram as this will interfere with the EM algorithm.
        for (unsigned int k = 0; k < smallSupport; k++)
                errorHist.erase(k);

        double edgeErrorLambda = edgeCovModel.getErrLambda();
        double edgeErrorODF = edgeCovModel.getErrorODF();
        double edgeErrorWeight = edgeCovModel.getWeight(0);

        // The below EM procedure might not convergence, but the obtained NB
        // parameters should be sensible in all cases: only small support
        // values are inferred, therefore, the mean should never be larger than
        // the mean of the error histogram (and thus be close to zero).
        nIter = Util::fitTruncNegBinomEM(errorHist, edgeErrorLambda,
                                         edgeErrorODF, edgeErrorWeight,
                                         epsilon, maxIter);

        if (nIter > maxIter)
                cout << "\tWARNING: EM algorithm to fit edge error model "
                        "did not converge\n";
        else
                cout << "\tEM algorithm to fit edge error model converged "
                     << "after " << nIter << " iteration(s)" << endl;

        // Compute edge average and weights
        vector<double> wEdge(edgeCovModel.getK(), 0.0);
        wEdge[0] = edgeErrorWeight;

        vector<double> totWEdgeCov(edgeCovModel.getK(),0.0);
        for (const auto& it : edgeMult) {
                Arc& arc = dBG.getArc(it.first);
                for (size_t m = 1; m < edgeCovModel.getK(); m++) {
                        wEdge[m] += it.second[m];
                        totWEdgeCov[m] += it.second[m] * arc.getCov();
                }
        }

        noml = 0.0, denl = 0.0;
        for (size_t m = 1; m < nodeCovModel.getK(); m++) {
                noml += totWEdgeCov[m] / m;
                denl += wEdge[m];
        }

        double edgeLambda = max(0.5, noml / denl);

        // Compute edge variance and ODF
        vector<double> wEdgeSqr(edgeCovModel.getK(), 0.0);

        for (const auto& it : edgeMult) {
                Arc& arc = dBG.getArc(it.first);
                for (size_t m = 1; m < edgeCovModel.getK(); m++) {
                        double delta = arc.getCov() - m * edgeLambda;
                        wEdgeSqr[m] += it.second[m] * delta * delta;
                }
        }

        nom = 0.0, den = 0.0;
        for (size_t m = 1; m < edgeCovModel.getK(); m++) {
                nom += wEdgeSqr[m] / (m * edgeLambda);
                den += wEdge[m];
        }
        double edgeODF = nom / den;
        edgeODF = max(1.0, edgeODF);

        edgeCovModel = CovModel(edgeErrorLambda, edgeErrorODF,
                                edgeLambda, edgeODF, wEdge);

        cout << "\tEdge spectrum: " << edgeCovModel << endl;

        return true;
}

int CRFMult::computeMultEM(NodeMap<Multiplicity>& nodeMult,
                           CovModel& nodeCovModel,
                           EdgeMap<Multiplicity>& edgeMult,
                           CovModel& edgeCovModel,
                           bool trainNodeODF,
                           double epsilon, int maxIter)
{
        int iter;
        for (iter = 1; iter <= maxIter; iter++) {
                cout << "Iteration " << iter << endl;
                // infer multiplicities given the model
                if (!Estep(nodeMult, nodeCovModel,
                           edgeMult, edgeCovModel, epsilon))
                                return iter;
                // update the model given the multiplicities
                if (!Mstep(nodeMult, nodeCovModel,
                           edgeMult, edgeCovModel, trainNodeODF, epsilon))
                                return iter;
        }

        return (iter);
}
