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

#include <sstream>
#include <random>
#include <utility>

#include "coverage.h"
#include "dbgraph.h"
#include "settings.h"
#include "util.h"
#include "alignment.h"
#include "nodeendstable.h"
#include "kmeroverlap.h"

using namespace std;

DBGraph::DBGraph(const DBGraph& dBG,
                 const NodeMap<int>& subNodes, const EdgeMap<int>& subEdges,
                 BiNodeMap<NodeID>& subToOrig, BiNodeMap<NodeID>& origToSub) :
        settings(dBG.settings), nodes(NULL), arcs(NULL), numNodes(0), numArcs(0)
{
        // +1 because index 0 is not used
        nodes = new DSNode[subNodes.size() + 1];

        // +2 because index 0 is not used, final index denotes 'end'.
        arcs = new Arc[2 * subEdges.size() + 2];

        // initialize the nodes
        subToOrig.clear();
        origToSub.clear();
        for (const auto [origNodeID, cov] : subNodes) {
                DSNode& subNode = getDSNode(++numNodes);

                subNode.setSequence(dBG.getDSNode(origNodeID).getSequence());
                subNode.setCov(subNode.getMarginalLength() * cov);

                subToOrig[ numNodes] =  origNodeID;
                subToOrig[-numNodes] = -origNodeID;
                origToSub[ origNodeID] =  numNodes;
                origToSub[-origNodeID] = -numNodes;
        }

        // initialize the arcs
        for (const auto [origNodeID, cov] : subNodes) {
                NodeID subNodeID = origToSub[origNodeID];
                DSNode& s = getDSNode(subNodeID);
                DSNode& o = dBG.getDSNode(origNodeID);

                s.setFirstLeftArc(arcs + numArcs + 1);
                for (ArcIt it = o.leftBegin(); it != o.leftEnd(); it++) {
                        EdgeRep er(it->getNodeID(), origNodeID);
                        if (subEdges.find(er) == subEdges.end())
                                continue;
                        NodeID leftID = origToSub[it->getNodeID()];
                        arcs[++numArcs].setNodeID(leftID);
                        arcs[numArcs].setCov(get(subEdges, er));
                        s.setNumLeftArcs(s.numLeftArcs() + 1);
                }

                s.setFirstRightArc(arcs + numArcs + 1);
                for (ArcIt it = o.rightBegin(); it != o.rightEnd(); it++) {
                        EdgeRep er(origNodeID, it->getNodeID());
                        if (subEdges.find(er) == subEdges.end())
                                continue;
                        NodeID rightID = origToSub[it->getNodeID()];
                        arcs[++numArcs].setNodeID(rightID);
                        arcs[numArcs].setCov(get(subEdges, er));
                        s.setNumRightArcs(s.numRightArcs() + 1);
                }
        }

        numValidNodes = numNodes;
        numValidArcs = numArcs;
}

void DBGraph::createNodeMapping(const DBGraph& sg, NodeID srcID, NodeID newSrcID,
                                BiNodeMap<NodePosPair>& mapping) const
{
        // given entry point for the mapping
        mapping[ srcID] = NodePosPair( newSrcID, 0);
        mapping[-srcID] = NodePosPair(-newSrcID, 0);

        vector<NodeID> stack({srcID, -srcID});
        BiNodeMap<bool> handled;

        while (!stack.empty()) {
                NodeID id = stack.back();
                stack.pop_back();

                if (handled[id])
                        continue;       // node already handled
                handled[id] = true;

                SSNode n = getSSNode(id);
                NodeLength len = n.getMarginalLength();

                auto [id_p, pos_p] = mapping[id];
                SSNode n_p = sg.getSSNode(id_p);
                NodeLength len_p = n_p.getMarginalLength();

                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        NodeID rID = it->getNodeID();
                        if (handled[rID])
                                continue;       // node already handled
                        SSNode r = getSSNode(rID);
                        NodeLength len_r = r.getMarginalLength();
                        char c = r.peekNucleotideMarginalLeft();

                        NodeID rID_p = 0; NodeLength pos_r_p = 0;
                        if (pos_p + len_p == len) {     // visit new node
                                rID_p = n_p.getRightArcNodeID(c);
                                if (rID_p == 0)
                                        continue;       // no such node
                        } else {                        // still in same node
                                if (n_p.getNucleotide(pos_p + len + Kmer::getK() - 1) != c)
                                        continue;
                                rID_p = id_p;
                                pos_r_p = pos_p + len;
                        }

                        SSNode r_p = sg.getSSNode(rID_p);
                        NodeLength len_r_p = r_p.getMarginalLength();
                        NodePosPair npp(rID_p, pos_r_p);
                        mapping[ rID] = npp;
                        npp = NodePosPair(-rID_p, len_r_p - pos_r_p - len_r);
                        mapping[-rID] = npp;

                        stack.push_back( rID);
                        stack.push_back(-rID);
                }
        }
}

bool DBGraph::consecutiveNPP(const NodePosPair& left,
                             const NodePosPair& right,
                             size_t offset) const
{
        // return false if one of the npps is invalid
        if (!left.isValid() || !right.isValid())
                return false;

        // left and right belong to the same node?
        if (left.getNodeID() == right.getNodeID())
                if ((left.getPosition() + offset) == right.getPosition())
                        return true;

        // left and right belong to connected nodes?
        SSNode lNode = getSSNode(left.getNodeID());
        if (lNode.rightArc(right.getNodeID()) == NULL)
                return false;

        // make sure the distance is consistent
        size_t dist = lNode.getMarginalLength() - left.getPosition() +
                                right.getPosition();
        return (offset == dist);
}

void DBGraph::writeBCalm (const string& filename)
{
        ofstream ofs(filename.c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename
                     << " for writing" << endl;

        for (NodeID id = 1; id <= numNodes; id++)
        {
                SSNode node = getSSNode(id);
                if (! node.isValid())
                        continue;

                ofs << ">" << id << " LN:i:" <<
                node.getMarginalLength() + Kmer::getK() -1 << " KC:f:" <<
                node.getCov() << " km:f:" << node.getAvgCov();
                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                        NodeID rightID = it->getNodeID();
                        ofs << "\tL:+:" << abs(rightID) << ":" << (rightID > 0 ? "+" : "-");
                }
                for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it ++) {
                        NodeID leftID = it->getNodeID();
                        ofs << "\tL:-:" << abs(leftID) << ":" << (leftID > 0 ? "-" : "+");
                }
                ofs << "\n";
                ofs << node.getSequence() << endl;
        }
        ofs.close();
}

void DBGraph::loadBCalm (const string& filename)
{
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw ios_base::failure("Cannot open file " + filename);

        // figure out file size
        ifs.seekg(0, ios_base::end);
        size_t fileSize = ifs.tellg();
        ifs.clear();
        ifs.seekg(0, ios::beg);

        // first pass through the file to find out number of nodes and arcs
        Util::startChrono();
        string progressStr = "Reading file " + filename + " (first pass)";

        numNodes = numArcs = 0;

        string line;
        while (getline(ifs, line)) {
                if (numNodes % 1024 == 0)
                        Util::progress(progressStr, ifs.tellg(), fileSize);
                if (line.front() != '>')
                        continue;

                numNodes++;

                // every occurrence of "L:" denotes a new arc
                size_t pos = 0;
                while ((pos = line.find("L:", pos + 1, 2)) != string::npos)
                        numArcs++;
        }

        double elapsed = Util::stopChrono();
        Util::progressEnd(progressStr, elapsed);

        // +1 because index 0 is not used
        nodes = new DSNode[numNodes+1];

        // +2 because index 0 is not used, final index denotes 'end'.
        arcs = new Arc[numArcs+2];

        // second pass through the file to store nodes and arcs
        cout << "Reading " << numNodes << " nodes and "
             << numArcs << " arcs..." << endl;

        ifs.clear();
        ifs.seekg(0, ios::beg);

        Util::startChrono();
        progressStr = "Reading file " + filename + " (second pass)";

        ArcID arcOffset = 1; NodeID nodeOffset = 1;
        while (getline(ifs, line)) {
                if (nodeOffset % 1024 == 0)
                        Util::progress(progressStr, ifs.tellg(), fileSize);
                DSNode& node = getDSNode(nodeOffset);

                // the line is a sequence
                if (line.front() != '>') {
                        node.setSequence(line);
                        nodeOffset++;
                        continue;
                }

                // the line is a FASTA descriptor line
                istringstream iss(line);

                // find "KC:i:" to figure out the k-mer coverage
                size_t pos = 0; int kmerCov = 0;
                if ((pos = line.find("KC:i:")) != string::npos) {
                        iss.seekg(pos + 5);
                        iss >> kmerCov;

                }

                node.setCov(kmerCov);

                vector<NodeID> leftArcs;
                vector<NodeID> rightArcs;

                // every occurrence of "L:" denotes a new arc
                pos = 0;
                while ((pos = line.find("L:", pos + 1, 2)) != string::npos) {
                        iss.seekg(pos + 2);
                        char c, l, r;
                        int dstID;
                        iss >> l >> c >> dstID >> c >> r >> right;
                        dstID++;        // we number from 1, not 0

                        if (l == '+')
                                rightArcs.push_back(r == '+' ? dstID : -dstID);
                        else
                                leftArcs.push_back(r == '-' ? dstID : -dstID);
                }

                node.setNumLeftArcs(leftArcs.size());
                node.setFirstLeftArc(arcs + arcOffset);
                for (auto dstID : leftArcs)
                        arcs[arcOffset++].setNodeID(dstID);

                node.setNumRightArcs(rightArcs.size());
                node.setFirstRightArc(arcs + arcOffset);
                for (auto dstID : rightArcs)
                        arcs[arcOffset++].setNodeID(dstID);
        }

        elapsed = Util::stopChrono();
        Util::progressEnd(progressStr, elapsed);
        ifs.close();

        // figure out the value for k
        bool autoDetectK = false;

        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;

                SSNode node = getSSNode(id);

                if (node.numRightArcs() < 2)
                        continue;

                ArcIt it = node.rightBegin();
                SSNode nb1 = getSSNode(it->getNodeID());
                it++;
                SSNode nb2 = getSSNode(it->getNodeID());
                string s1 = nb1.getSequence();
                string s2 = nb2.getSequence();

                int k = 0;
                while (s1[k] == s2[k])          // overlap == k-1
                        k++;
                k++;                            // add 1 to get k

                Kmer::setWordSize(k);
                autoDetectK = true;
                break;
        }

        if (!autoDetectK)
                throw runtime_error("Cannot infer kmer size from input file");
        else
                cout << "Kmer size is " << Kmer::getK() << endl;

        numValidNodes = numNodes;
        numValidArcs = numArcs;
}

void DBGraph::writeBinary(const std::string& filename) const
{
        ofstream ofs(filename.c_str(), ios::binary);

        // A) write the header
        size_t k = Kmer::getK();
        ofs.write((char*)&k, sizeof(k));
        ofs.write((char*)&numNodes, sizeof(numNodes));
        ofs.write((char*)&numArcs, sizeof(numArcs));

        // B) write nodes
        for (NodeID id = 1; id <= numNodes; id++)
                getDSNode(id).write(ofs, arcs);

        // C) write arcs
        for (ArcID i = 0; i < numArcs+2; i++)
                arcs[i].write(ofs);

        ofs.close();

        cout << "Wrote " << numNodes << " nodes and "
             << numArcs << " arcs" << endl;
}

void DBGraph::loadBinary(const std::string& filename)
{
        // read the metadata
        ifstream ifs(filename.c_str(), ios::binary);
        if (!ifs)
                throw ios_base::failure("Cannot open " + filename);

        // A) read the header
        size_t k;
        ifs.read((char*)&k, sizeof(k));
        Kmer::setWordSize(k);
        ifs.read((char*)&numNodes, sizeof(numNodes));
        ifs.read((char*)&numArcs, sizeof(numArcs));
        // +1 because index 0 isn't used
        nodes = new DSNode[numNodes+1];
        // +2 because index 0 isn't used, final index denotes 'end'
        arcs = new Arc[numArcs+2];

        // B) create the nodes
        for (NodeID id = 1; id <= numNodes; id++)
                getDSNode(id).read(ifs, arcs);

        // C) create the arcs
        for (ArcID i = 0; i < numArcs+2; i++)
                arcs[i].read(ifs);

        ifs.close();

        numValidNodes = numNodes;
        numValidArcs = numArcs;
}

void DBGraph::loadGraph(const string& filename)
{
        ifstream ifsMeta(filename + ".meta.st1");
        if (!ifsMeta)
                throw ios_base::failure("Cannot open file " + filename);
        size_t k;
        ifsMeta >> k >> numNodes >> numArcs;
        ifsMeta.close();

        Kmer::setWordSize(k);

        // +1 because index 0 is not used
        nodes = new DSNode[numNodes+1];

        // +2 because index 0 is not used, final index denotes 'end'.
        arcs = new Arc[numArcs+2];

        cout << "Reading de Bruijn graph (k = " << Kmer::getK() << ") with "
             << numNodes << " nodes and " << numArcs << " edges..." << endl;

        // read the nodes
        Util::startChrono();
        string progressStr = "Reading file " + filename + ".node.st1";
        ifstream ifsNode(filename + ".node.st1");
        if (!ifsNode)
                throw ios_base::failure("Cannot open file " + filename);

        NodeEndTable table(true, 2 * numNodes);

        for (size_t i = 1; i <= numNodes; i++) {
                if (i % 1024 == 0)
                        Util::progress(progressStr, i, numNodes);

                string seq;
                ifsNode >> seq;
                if (!ifsNode)
                        break;

                DSNode& node = getDSNode(i);
                node.setSequence(seq);

                Kmer firstKmer = node.getLeftKmer();
                Kmer finalKmer = node.getRightKmer();

                if (node.getMarginalLength() > 1) {
                        if (!table.insert(firstKmer, i))
                                cerr << "ERROR: Multiple nodes start/end with the same k-mer!" << endl;
                        if (!table.insert(finalKmer, i))
                                cerr << "ERROR: Multiple nodes start/end with the same k-mer!" << endl;
                } else if (!table.insert(firstKmer, i))
                        cerr << "ERROR: Multiple nodes start/end with the same k-mer!" << endl;
        }

        ifsNode.close();
        double elapsed = Util::stopChrono();
        Util::progressEnd(progressStr, elapsed);

        // read the arcs
        Util::startChrono();
        progressStr = "Reading file " + filename + ".edge.st1";
        ifstream ifsEdge(filename + ".edge.st1");
        if (!ifsEdge)
                throw ios_base::failure("Cannot open file " + filename);

        size_t arcOffset = 1;
        for (size_t i = 1; i <= numNodes; i++) {
                if (i % 1024 == 0)
                        Util::progress(progressStr, i, numNodes);

                int bfLeft, bfRight;
                ifsEdge >> bfLeft >> bfRight;
                if (!ifsEdge)
                        break;

                // read the arc info
                KmerOverlap overlap((bfLeft << 4) + bfRight);

                int numLeftArcs = overlap.getNumLeftOverlap();
                int numRightArcs = overlap.getNumRightOverlap();

                DSNode& node = getDSNode(i);

                node.setNumLeftArcs(numLeftArcs);
                node.setNumRightArcs(numRightArcs);

                node.setFirstLeftArc(arcs + arcOffset);
                node.setFirstRightArc(arcs + arcOffset + numLeftArcs);

                // connect the left arcs alphabetically (ACGT)
                Kmer firstKmer = node.getLeftKmer();
                for (NucleotideID j = 0; j < 4; j++) {
                        char n = Nucleotide::nucleotideToChar(j);
                        if (!overlap.hasLeftOverlap(n))
                                continue;

                        Kmer kmer = firstKmer;
                        kmer.pushNucleotideLeft(n);

                        NodeEndRef ref = table.find(kmer);
                        if (ref.first == table.end())
                                throw ios_base::failure("Mismatch between nodes"
                                "and arc file.");
                        arcs[arcOffset++].setNodeID(ref.getNodeID());
                }

                // connect the right arcs alphabetically (ACGT)
                Kmer finalKmer = node.getRightKmer();
                for (NucleotideID j = 0; j < 4; j++) {
                        char n = Nucleotide::nucleotideToChar(j);
                        if (!overlap.hasRightOverlap(n))
                                continue;

                        Kmer kmer = finalKmer;
                        kmer.pushNucleotideRight(n);

                        NodeEndRef ref = table.find(kmer);
                        if (ref.first == table.end())
                                throw ios_base::failure("Mismatch between nodes"
                                "and arc file.");
                        arcs[arcOffset++].setNodeID(ref.getNodeID());
                }
        }

        ifsEdge.close();
        elapsed = Util::stopChrono();
        Util::progressEnd(progressStr, elapsed);

        numValidNodes = numNodes;
        numValidArcs = numArcs;
}

void DBGraph::writeContigs(const string& filename,
                           const vector<NodeChain>& contigs,
                           const vector<string>& contigNames) const
{
        assert(contigs.size() == contigNames.size());

        ofstream ofs(filename.c_str());
        for (size_t i = 0; i < contigs.size(); i++) {
                if (contigs[i].empty())
                        continue;
                string str = convertNodesToString(contigs[i]);

                ofs << ">" << contigNames[i] << "\n";
                Util::writeSeqWrap(ofs, str, 60);
        }
}

void DBGraph::getGraph(std::vector<NodeID>& nodes, std::vector<EdgeID>& edges)
{
        nodes.reserve(2 * getNumValidNodes());
        edges.reserve(getNumValidArcs());

        for (NodeID srcID = -numNodes; srcID <= numNodes; srcID++) {
                if (srcID == 0)
                        continue;
                SSNode n = getSSNode(srcID);
                if (!getSSNode(srcID).isValid())
                        continue;
                nodes.push_back(srcID);

                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        NodeID dstID = it->getNodeID();
                        edges.push_back(make_pair(srcID, dstID));
                }
        }
}

void DBGraph::getSubgraphAux(
        priority_queue<NodeDepth, vector<NodeDepth>, NodeDepthComp>& todo,
        vector<NodeID>& nodes, vector<EdgeID>& edges, size_t maxDepth) const
{
        set<NodeID> nodeSet;

        while (!todo.empty()) {
                // get and erase the current node
                auto [thisID, thisDepth] = todo.top();
                todo.pop();

                if (!nodeExists(thisID))
                        continue;

                SSNode n = getSSNode(thisID);

                // if the node was already handled, skip
                if (nodeSet.find(thisID) != nodeSet.end())
                        continue;

                nodeSet.insert(thisID);

                // don't go any deeper if we've reached the maximum depth
                if (thisDepth >= maxDepth)
                        continue;

                // process the right arcs
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        NodeID rightID = it->getNodeID();
                        if (nodeSet.find(rightID) != nodeSet.end())
                                continue;       // edge already added by rightID

                        edges.push_back(make_pair(thisID, rightID));
                        todo.push(NodeDepth(rightID, thisDepth+1));
                }

                // process the left arcs
                for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++) {
                        NodeID leftID = it->getNodeID();
                        if (nodeSet.find(leftID) != nodeSet.end())
                                continue;       // edge already added by leftID

                        edges.push_back(make_pair(leftID, thisID));
                        todo.push(NodeDepth(leftID, thisDepth + 1));
                }
        }

        nodes = vector<NodeID>(nodeSet.begin(), nodeSet.end());
}

void DBGraph::getSubgraph(NodeID seedNode, vector<NodeID>& nodes,
                          vector<EdgeID>& edges, size_t maxDepth) const
{
        // initialize priority queue with seed node
        priority_queue<NodeDepth, vector<NodeDepth>, NodeDepthComp> pq;
        pq.push(NodeDepth(seedNode, 0));

        // run a search
        getSubgraphAux(pq, nodes, edges, maxDepth);
}

void DBGraph::getSubgraph(const NodeChain& seedChain, vector<NodeID>& nodes,
                          vector<EdgeID>& edges, size_t maxDepth) const
{
        // initialize priority queue with nodes from the node chain
        priority_queue<NodeDepth, vector<NodeDepth>, NodeDepthComp> pq;
        for (NodeID id : seedChain)
                pq.push(NodeDepth(id, 0));

        // run a search
        getSubgraphAux(pq, nodes, edges, maxDepth);
}

vector<NodeRep> DBGraph::getNodeReps(const string& filename) const
{
        ifstream ifs(filename.c_str());

        NodeID n;
        size_t numNodes = 0;
        while (ifs >> n)
                numNodes++;

        vector<NodeRep> nodes;
        nodes.reserve(numNodes);

        ifs.clear();
        ifs.seekg(0, ios::beg);

        while (ifs >> n) {
                if (!nodeExists(n))
                        throw runtime_error("Node with ID " + to_string(n) + " does not exist");
                nodes.push_back(NodeRep(n));
        }

        return nodes;
}

vector<NodeRep> DBGraph::getNodeReps(size_t N) const
{
        vector<NodeRep> nodes;
        nodes.reserve(getNumValidNodes());
        for (size_t i = 1; i <= numNodes; i++)
                if (getSSNode(i).isValid())
                        nodes.push_back(NodeRep(i));

        N = min(N, nodes.size());
        if (N == nodes.size())  // if you need all nodes we do not shuffle
                return nodes;

        // sample N nodes using the Fisher-Yates algoritm
        random_device rd;
        mt19937 mt;     // uses the default, fixed seed for reproducibility
        for (size_t i = 0; i < N; i++) {
                uniform_int_distribution<size_t> dis(i, nodes.size() - 1);
                swap(nodes[i], nodes[dis(mt)]);
        }

        return vector<NodeRep>(nodes.begin(), nodes.begin() + N);
}

vector<EdgeRep> DBGraph::getEdgeReps(const string& filename) const
{
        ifstream ifs(filename.c_str());

        ArcID l, r;
        size_t numEdges = 0;
        while (ifs >> l >> r)
                numEdges++;

        vector<EdgeRep> edges;
        edges.reserve(numEdges);

        ifs.clear();
        ifs.seekg(0, ios::beg);

        while (ifs >> l >> r) {
                if (!edgeExists(EdgeRep(l, r)))
                        throw runtime_error("Edge from ID " + to_string(l) + " to ID " + to_string(r) + " does not exist");
                edges.push_back(EdgeRep(l, r));
        }

        return edges;
}

vector<EdgeRep> DBGraph::getEdgeReps(size_t N) const
{
        vector<EdgeRep> edges;
        edges.reserve(getNumValidArcs() / 2);   // roughly half of the nodes
        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;
                SSNode n = getSSNode(id);
                if (!n.isValid())
                        continue;
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++)
                        if (id <= -it->getNodeID())
                                edges.push_back(EdgeRep(id, it->getNodeID()));
        }

        N = min(N, edges.size());
        if (N == edges.size())  // if you need all edges we do not shuffle
                return edges;

        // sample N edges using the Fisher-Yates algoritm
        random_device rd;
        mt19937 mt;     // uses the default, fixed seed for reproducibility
        for (size_t i = 0; i < N; i++) {
                uniform_int_distribution<size_t> dis(i, edges.size() - 1);
                swap(edges[i], edges[dis(mt)]);
        }

        return vector<EdgeRep>(edges.begin(), edges.begin() + N);
}

vector<EdgeRep> DBGraph::getLowCovEdges(double threshold) const
{
        vector<EdgeRep> edges;
        edges.reserve(getNumValidArcs() / 2);   // roughly half of the nodes
        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;
                SSNode n = getSSNode(id);
                if (!n.isValid())
                        continue;
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        if (id <= -it->getNodeID() && it->getCov() <= threshold)
                                edges.push_back(EdgeRep(id, it->getNodeID()));
                }
        }

        return edges;
}

vector<NodeRep> DBGraph::getLowCovForks(double threshold) const
{
        /**
         * Fork has 1 left neighbor and >= 1 parallel nodes with higher coverage
         *         |--> p1
         * root -->|--> fork  <- node with lowest coverage
         *         |--> p2
         */
        vector<NodeRep> nodes;
        for (NodeID i = -numNodes; i <= numNodes; i++) {
                if (i == 0)
                        continue;
                SSNode n = getSSNode(i);
                if (!n.isValid())
                        continue;

                if (n.getAvgCov() > threshold)
                        continue;       // coverage too high

                if (n.numLeftArcs() != 1)
                        continue;       // not a fork

                SSNode l = getSSNode(n.leftBegin()->getNodeID());
                if (n.getAvgCov() > l.getAvgCov())
                        continue;       // coverage too high

                if (l.numRightArcs() < 2)
                        continue;       // not a fork

                /*if (l.numLeftArcs() != 1)
                        continue;       // not a fork*/

                bool isFork = true;
                for (ArcIt it = l.rightBegin(); it != l.rightEnd(); it++) {
                        if (it->getNodeID() == i)
                                continue;
                        SSNode p = getSSNode(it->getNodeID());
                        if (p.getAvgCov() > l.getAvgCov())
                                isFork = false; // not a fork
                        if (n.getAvgCov() > p.getAvgCov())
                                isFork = false; // coverage too high
                }

                if (isFork)
                        nodes.push_back(NodeRep(i));
        }

        return nodes;
}

std::vector<NodeRep> DBGraph::getLowCovForks(double loCov, double hiCov) const
{
        /**
         * Fork has 1 left neighbor and >= 1 parallel nodes with higher coverage
         *         |--> p1
         * root -->|--> fork  <- node with lowest coverage
         *         |--> p2
         */
        vector<NodeRep> nodes;
        for (NodeID i = -numNodes; i <= numNodes; i++) {
                if (i == 0)
                        continue;
                SSNode n = getSSNode(i);
                if (!n.isValid())
                        continue;

                if (n.getAvgCov() > loCov)
                        continue;       // coverage too high

                if (n.numLeftArcs() != 1)
                        continue;       // not a fork

                SSNode l = getSSNode(n.leftBegin()->getNodeID());
                if (n.getAvgCov() > l.getAvgCov())
                        continue;       // coverage too high

                /*if (l.numLeftArcs() != 1)
                        continue;       // not a fork*/

                if (l.numRightArcs() < 2)
                        continue;       // not a fork

                if (l.getAvgCov() < hiCov)
                        continue;       // coverage too low

                bool isFork = true;
                for (ArcIt it = l.rightBegin(); it != l.rightEnd(); it++) {
                        if (it->getNodeID() == i)
                                continue;
                        SSNode p = getSSNode(it->getNodeID());
                        if (p.getAvgCov() > l.getAvgCov())
                                isFork = false; // not a fork
                        if (n.getAvgCov() > p.getAvgCov())
                                isFork = false; // coverage too high
                }

                if (isFork)
                        nodes.push_back(NodeRep(i));
        }

        return nodes;
}

std::vector<NodeRep> DBGraph::getNonInterfering(double threshold)
{
        vector<pair<double, NodeRep>> cnp;
        for (size_t i = 1; i <= numNodes; i++) {
                SSNode n = getSSNode(i);
                if (!n.isValid())
                        continue;

                if (n.getAvgCov() <= threshold)
                        cnp.push_back(make_pair(n.getAvgCov(), NodeRep(i)));
        }

        cout << "Total number of candidate nodes: " << cnp.size() << endl;
        sort(cnp.begin(), cnp.end());

        vector<NodeRep> nodes;
        for (const auto& el : cnp) {
                NodeRep nr = el.second;
                if (getSSNode(nr).getFlag2())
                        continue;
                vector<NodeID> sgn; vector<EdgeID> sge;
                getSubgraph(nr, sgn, sge, 3);

                bool graphFlagged = false;
                for (const auto& id : sgn) {
                        if (getSSNode(id).getFlag1()) {
                                graphFlagged = true;
                                break;
                        }
                }

                if (graphFlagged)
                        continue;

                nodes.push_back(nr);
                for (const auto& id : sgn)
                        getSSNode(id).setFlag1(true);
        }

        setAllFlags1(false);
        setAllFlags2(false);

        return nodes;
}

std::vector<NodeRep> DBGraph::getLowCovNodesAdjToHighCov(double loCov,
                                                         double hiCov) const
{
        vector<NodeRep> nodes;
        nodes.reserve(getNumValidNodes());
        for (size_t i = 1; i <= numNodes; i++) {
                SSNode n = getSSNode(i);
                if (!n.isValid())
                        continue;

                if (n.getAvgCov() > loCov)
                        continue;

                if ((n.numLeftArcs() != 1) || (n.numRightArcs() != 1))
                        continue;       // not a bubble

                if (n.getMarginalLength() > 2*Kmer::getK())
                        continue;       // not a bubble

                if ((getSSNode(n.leftBegin()->getNodeID()).getAvgCov() < hiCov) ||
                    (getSSNode(n.rightBegin()->getNodeID()).getAvgCov() < hiCov))
                        continue;       // coverage too high

                nodes.push_back(NodeRep(i));
        }

        return nodes;
}

vector<NodeID> DBGraph::findParallelPath(NodeID id)
{
        vector<NodeID> result;
        SSNode noi = getSSNode(id);     // node-of-interest

        // dstID = most downstream, non-branching node
        NodeID dstID = noi.rightBegin()->getNodeID();
        size_t rlen = 0;
        while (true) {
                SSNode dst = getSSNode(dstID);
                if (dst.numRightArcs() != 1)
                        break;
                if (rlen + dst.getMarginalLength() > MAX_BUBBLE_LEN)
                        break;
                dstID = dst.rightBegin()->getNodeID();
                rlen += dst.getMarginalLength();
        }

        priority_queue<NodeDepth, vector<NodeDepth>, NodeDepthComp> pq;
        BiNodeMap<NodeID> prev;
        BiNodeMap<int> dist;
        set<NodeID> visited;

        // we consider multiple srcID, i.e., all upstream non-branching nodes
        NodeID srcID = noi.leftBegin()->getNodeID();
        size_t llen = 0;
        do {
                pq.push(NodeDepth(srcID, llen));        // add to pq
                prev[srcID] = 0;
                dist[srcID] = llen;

                SSNode src = getSSNode(srcID);
                if (src.numLeftArcs() != 1)
                        break;
                srcID = src.leftBegin()->getNodeID();
                llen += src.getMarginalLength();
        } while (llen <= MAX_BUBBLE_LEN);

        size_t maxLen = rlen + 2*llen + noi.getMarginalLength() + 3;
        int counter = 0;

        while (!pq.empty())
        {
                // get element with smallest currDist from the priority queue
                auto [currID, currDist] = pq.top();
                pq.pop();

                if (!visited.insert(currID).second)     // mark as visited
                        continue;       // node already visited

                if (currID == dstID) {  // extract path and exit
                        while (currID != 0) {
                                result.push_back(currID);
                                currID = prev[currID];
                        }
                        reverse(result.begin(), result.end());
                        return result;
                }

                if (counter++ > 100)
                        return vector<NodeID>();

                if (currDist > maxLen)
                        continue;       // maximum length exceeded

                SSNode n = getSSNode(currID);
                if (n.getAvgCov() <= noi.getAvgCov())
                        continue;       // insufficient coverage

                // process the right arcs
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        NodeID rID = it->getNodeID();
                        int d = currDist + getSSNode(rID).getMarginalLength();
                        pq.push(NodeDepth(rID, d));
                        if ((dist.find(rID) == dist.end()) || (d < dist[rID])) {
                                dist[rID] = d;
                                prev[rID] = currID;
                        }
                }
        }

        /*if (visited.size() > 10000) {
                cout << visited.size() << endl;
                cout << maxLen << endl;
        }*/

        return result;
}

vector<pair<NodeID, double>> DBGraph::getRightNeigbors(NodeID id) const
{
        SSNode n = getSSNode(id);

        vector<pair<NodeID, double>> result;
        result.reserve(n.numRightArcs());

        for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                NodeID rID = it->getNodeID();
                result.emplace_back(rID, getSSNode(rID).getAvgCov());
        }

        // sort the next nodes by coverage (low to high)
        auto sortByCov = [] (pair<NodeID, double> const& s1,
                             pair<NodeID, double> const& s2) -> bool {
                return s1.second < s2.second;
        };
        sort(result.begin(), result.end(), sortByCov);

        return result;
}

void DBGraph::smoothBubbles(size_t maxLen, size_t maxED)
{
        vector<NodeRep> nodes;
        for (size_t i = 1; i <= numNodes; i++) {
                SSNode n = getSSNode(i);
                if (!n.isValid())
                        continue;       // deleted node

                if ((n.numLeftArcs() != 1) || (n.numRightArcs() != 1))
                        continue;       // not a bubble

                if (n.getMarginalLength() > maxLen)
                        continue;       // not a bubble

                // an alternative path should exist between leftID and rightID
                // where each node in the path has avgCov > n.avgCov.
                NodeID leftID = n.leftBegin()->getNodeID();
                NodeID rightID = n.rightBegin()->getNodeID();

                if (abs(leftID) == abs(rightID))
                        continue;       // skip loops and palindromes

                if (!findPath(leftID, rightID, n.getAvgCov(), maxLen))
                        continue;

                nodes.push_back(NodeRep(i));
        }

        cout << "Removing " << nodes.size() << " nodes\n";
        removeNodes(nodes);
}

vector<NodeRep> DBGraph::getLowCovChimConn(double threshold,
                                           size_t maxLen) const
{
        vector<NodeRep> nodes;
        for (size_t i = 1; i <= numNodes; i++) {
                SSNode n = getSSNode(i);
                if (!n.isValid())
                        continue;       // deleted node

                if ((n.numLeftArcs() != 1) || (n.numRightArcs() != 1))
                        continue;       // not a bubble

                if (n.getMarginalLength() > maxLen)
                        continue;       // node too long

                if (n.getAvgCov() > threshold)
                        continue;       // coverage too high

                if (getSSNode(n.leftBegin()->getNodeID()).getAvgCov() < 4.0 * n.getAvgCov())
                        continue;

                if (getSSNode(n.rightBegin()->getNodeID()).getAvgCov() < 4.0 * n.getAvgCov())
                        continue;

                nodes.push_back(NodeRep(i));
        }

        return nodes;
}

void DBGraph::getNodeCovHist(map<int, double>& hist)
{
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;

                hist[node.getAvgCov() + 0.5]++;// node.getMarginalLength();
        }
}

void DBGraph::writeCytoscapeGraph(const std::string& filename,
                                  vector<NodeID> nodes,
                                  vector<EdgeID> edges,
                                  const NodeMap<Multiplicity>& enm,
                                  const EdgeMap<Multiplicity>& eem,
                                  const NodeMap<int>& tnm,
                                  const EdgeMap<int>& tem) const
{
        // A) write all arcs
        ofstream ofs((filename + ".arcs").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".arcs"
                     << " for writing" << endl;
        ofs << "Source node\tTarget node\tCoverage\tTrue mult.\t"
               "Est mult.\tLOR\tBelief\tcorrect\n";

        for (const auto& edge : edges) {
                EdgeRep er(edge);
                ArcID id = getArcID(er);

                auto teIt = tem.find(er);
                int trueMult = (teIt == tem.end()) ? -1 : teIt->second;

                auto eeIt = eem.find(er);
                int estMult = (eeIt == eem.end()) ? -1 : eeIt->second.getExpMult();

                double multLOR = (eeIt == eem.end()) ?
                        -1.0 : eeIt->second.getExpMultLogOR();
                double multBelief = (eeIt == eem.end()) ?
                        -1.0 : exp(eeIt->second.getExpMultLProb());

                ofs << edge.first << "\t" << edge.second << "\t"
                    << getArc(id).getCov() << "\t" << trueMult << "\t"
                    << estMult << "\t" << multLOR << "\t" << multBelief << "\t" << (trueMult==estMult? 1:0)<< "\n";
        }
        ofs.close();

        // B) write all nodes
        ofs.open((filename + ".nodes").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".nodes"
                     << " for writing" << endl;

        ofs << "Node ID\tMarginal length\tLeft arcs\tRight arcs\tCoverage\t"
               "True mult.\tEst. mult.\tLOR\tBelief\tcorrect\tSequence\n";

        for (const auto& id : nodes) {
                NodeRep nr(id);
                SSNode n = getSSNode(id);

                auto tnIt = tnm.find(nr);
                int trueMult = (tnIt == tnm.end()) ? -1 : tnIt->second;

                auto enIt = enm.find(nr);
                int estMult = (enIt == enm.end()) ? -1 : enIt->second.getExpMult();

                double multLOR = (enIt == enm.end()) ?
                        -1.0 : enIt->second.getExpMultLogOR();
                double multBelief = (enIt == enm.end()) ?
                        -1.0 : exp(enIt->second.getExpMultLProb());

                ofs << id << "\t" << n.getMarginalLength() << "\t"
                    << (int)n.numLeftArcs() << "\t" << (int)n.numRightArcs()
                    << "\t" << n.getAvgCov() << "\t" << trueMult << "\t"
                    << estMult << "\t" << multLOR << "\t" << multBelief<< "\t" << (trueMult==estMult? 1:0) << "\t"
                    << n.getSequence() << "\n";
        }
        ofs.close();
}

void DBGraph::writeCytoscapeSubgraph(const std::string& filename,
                                     vector<NodeID> nodes,
                                     vector<EdgeID> edges,
                                     BiNodeMap<NodeID> subToOrig) const
{
        // A) write all arcs
        ofstream ofs((filename + ".arcs").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".arcs"
                     << " for writing" << endl;
        ofs << "Source node\tTarget node\tCoverage\n";

        for (const auto& edge : edges) {
                EdgeRep er(edge);
                ArcID id = getArcID(er);

                ofs << subToOrig[edge.first] << "\t"
                   << subToOrig[edge.second] << "\t"
                    << getArc(id).getCov() << "\n";
        }
        ofs.close();

        // B) write all nodes
        ofs.open((filename + ".nodes").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".nodes"
                     << " for writing" << endl;

        ofs << "Node ID\tMarginal length\tLeft arcs\tRight arcs\tCoverage\t"
               "Sequence\n";

        for (const auto& id : nodes) {
                NodeRep nr(id);
                SSNode n = getSSNode(id);

                ofs << subToOrig[id] << "\t" << n.getMarginalLength() << "\t"
                    << (int)n.numLeftArcs() << "\t" << (int)n.numRightArcs()
                    << "\t" << n.getAvgCov() << "\t" << n.getSequence() << "\n";
        }
        ofs.close();
}

void DBGraph::defragNodes()
{
        vector<NodeID> old2new(numNodes + 1, 0);

        // defragment the node array
        for (NodeID oldID = 1, newID = 1; oldID <= numNodes; oldID++) {
                SSNode n = getSSNode(oldID);
                if (!n.isValid())
                        continue;

                old2new[oldID] = newID;

                // we use the move assignment operator for efficiency
                // (to avoid a deep copy of the TString member)
                nodes[newID] = move(nodes[oldID]);
                newID++;
        }

        // update the arcs to point to the new nodeIDs
        for (ArcID id = 1; id <= numArcs; id++) {
                NodeID oldID = arcs[id].getNodeID();
                NodeID newID = (oldID < 0) ? -old2new[-oldID] : old2new[oldID];
                arcs[id].setNodeID(newID);
        }

        numNodes = numValidNodes;
        cout << "Defragged nodes array: " << numNodes << " nodes\n";
}

void DBGraph::defragArcs()
{
        vector<ArcID> old2new(numArcs + 2, 0);

        numValidArcs = 0;
        for (ArcID oldID = 1; oldID <= numArcs; oldID++) {
                if (!arcs[oldID].isValid())
                        continue;

                numValidArcs++;
                old2new[oldID] = numValidArcs;

                // we use the regular assignment operator
                arcs[numValidArcs] = arcs[oldID];
        }

        // update the nodes to point to the new arcIDs
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode n = getSSNode(id);
                if (!n.isValid())
                        continue;

                ArcID oldID = n.getFirstLeftArc() - arcs;
                ArcID newID = old2new[oldID];
                n.setFirstLeftArc(arcs + newID);

                oldID = n.getFirstRightArc() - arcs;
                newID = old2new[oldID];
                n.setFirstRightArc(arcs + newID);
        }

        numArcs = numValidArcs;
        cout << "Defragged arcs array: " << numArcs << " arcs\n";
}

void DBGraph::createKmerNPPTable(KmerNPPTable& table) const
{
        //Util::startChrono();
        //string progressStr("Populating <k-mer, node> table");

        // count the number of k-mers in the graph
        size_t numKmers = 0;
        for (NodeID id = 1; id <= numNodes; id++) {
                const DSNode& node = nodes[id];
                if (!node.isValid())
                        continue;
                numKmers += node.getMarginalLength();
        }

        table.resize(numKmers);

        // populate the table with kmers
        for (NodeID id = 1; id <= numNodes; id++) {
                /*if (id % 1024 == 0)
                        Util::progress(progressStr, id, numNodes);*/

                const DSNode &node = nodes[id];
                if (!node.isValid())
                        continue;
                const TString& tStr = node.getTSequence();
                Kmer kmer(tStr);
                NodePosition pos = 0;
                table.insert(kmer, NodePosPair(id, pos++));

                for (size_t i = Kmer::getK(); i < tStr.getLength(); i++) {
                        kmer.pushNucleotideRight(tStr[i]);
                        table.insert(kmer, NodePosPair(id, pos++));
                }
        }

        //double elapsed = Util::stopChrono();
        //Util::progressEnd(progressStr, elapsed);
}

int DBGraph::getAbundanceMin() const {
        return settings.getAbundanceMin();
}

void DBGraph::getContigs(const LibraryContainer& libraries,
                         NodeMap<int>& nodeMult, EdgeMap<int>& edgeMult,
                         vector<NodeChain>& contigs,
                         vector<string>& contigName)
{

        vector<NodeChain> reductions, looseEnds;
        getReductions(libraries, nodeMult, edgeMult, 0.25f, 2,
                      reductions, looseEnds);

        vector<NodeChain> filteredLE;
        for (const auto& e : looseEnds) {
                bool keepLE = true;
                for (size_t i = 1; i < e.size(); i++)
                        if (nodeMult[e[i]] < 2)
                                keepLE = false;
                if (keepLE)
                        filteredLE.push_back(e);
        }
        looseEnds = filteredLE;
        filteredLE.clear();

        /*for (auto e : looseEnds)
                cout << e << endl;*/

        sort(looseEnds.begin(), looseEnds.end(), NodeChain::covHiToLo);

        BiNodeMap<bool> usedRepeats; BiNodeMap<NodeChain> filtLE;
        for (const auto& nc : looseEnds) {
                // we do not want big repeat nodes at both sides of its
                // flanking unique contigs
                if (usedRepeats[-nc.back()] == true)
                        continue;
                for (size_t i = 1; i < nc.size(); i++)
                        usedRepeats[nc[i]] = true;
                filtLE[nc.front()] = nc;
        }

        for (const auto& e : filtLE) {
                NodeID id = e.first;
                NodeChain m;
                if (filtLE.find(-id) != filtLE.end()) {
                        if (id > -id)
                                continue;       // already handled
                        m = NodeChain::merge(filtLE[-id].getRevCompl(), e.second);
                } else {
                        m = e.second;
                }
                filteredLE.push_back(m);
        }

        vector<pair<size_t, NodeChain>> lenID;

        // convert loose ends to contigs
        setAllFlags1(false);
        for (const auto& nc : filteredLE) {
                if (nc.empty())
                        continue;
                for (NodeID id : nc)
                        getSSNode(id).setFlag1(true);

                size_t len = 0;
                for (NodeID id : nc)
                        len += getSSNode(id).getMarginalLength();
                lenID.push_back(make_pair(len, nc));
        }

        // add the other nodes
        for (NodeID id = 1; id <= getNumNodes(); id++) {
                SSNode n = getSSNode(id);
                if (!n.isValid() || n.getFlag1())
                        continue;
                lenID.push_back(make_pair(n.getMarginalLength(),
                        NodeChain(vector<NodeID>(1, id))));
        }
        setAllFlags1(false);

        // sort large-to-small (C++14)
        sort(lenID.begin(), lenID.end(), std::greater<>());

        // produce contigs and contig names
        contigs.resize(lenID.size());
        contigName.resize(lenID.size());

        for (size_t i = 0; i < lenID.size(); i++) {
                NodeID id = lenID[i].second.front();
                SSNode n = getSSNode(id);
                contigs[i] = lenID[i].second;
                contigName[i] = "CONTIG_" + to_string(i+1) +
                                "_node_" + to_string(id) +
                                "_len_" + to_string(lenID[i].first) +
                                "_cov_" + to_string(n.getAvgCov());
        }
}

void DBGraph::getSingletonContigs(NodeMap<Multiplicity>& nodeMult,
                                  vector<NodeChain>& contigs,
                                  vector<string>& contigName) const
{
        // temporary vector to sort nodes by length
        vector<pair<size_t, NodeID>> lenID;
        lenID.reserve(numNodes);
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode n = getSSNode(id);
                if (!n.isValid())
                        continue;
                lenID.push_back(make_pair(n.getMarginalLength(), id));
        }

        // sort large-to-small (C++14)
        sort(lenID.begin(), lenID.end(), std::greater<>());

        // produce contigs and contig names
        contigs.reserve(lenID.size());
        contigName.reserve(lenID.size());

        size_t contigID = 1;
        for (auto [ignore, id] : lenID) {
                SSNode n = getSSNode(id);
                size_t M = nodeMult[id].getExpMult();
                contigs.emplace_back(NodeChain(vector<NodeID>(1, id)));
                contigName.emplace_back("CONTIG_" + to_string(contigID++) +
                                        "_node_" + to_string(id) +
                                        "_len_" + to_string(n.length()) +
                                        "_cov_" + to_string(n.getAvgCov()) +
                                        "_mult_" + to_string(M));
        }
}

NodeChain DBGraph::getConcatenation(NodeID seedID)
{
        NodeChain result;

        SSNode seed = getSSNode(seedID);
        if (!seed.isValid() || seed.getFlag1())
                return result;

        deque<NodeID> nodeListq;
        nodeListq.push_back(seedID);
        seed.setFlag1(true);

        // find linear paths to the right
        SSNode curr = seed;
        while (curr.numRightArcs() == 1) {
                NodeID rightID = curr.rightBegin()->getNodeID();
                SSNode right = getSSNode(rightID);
                // don't merge palindromic repeats / loops
                if (right.getFlag1())
                        break;
                if (right.numLeftArcs() != 1)
                        break;
                nodeListq.push_back(rightID);
                right.setFlag1(true);
                curr = right;
        }

        // find linear paths to the left
        curr = seed;
        while (curr.numLeftArcs() == 1) {
                NodeID leftID = curr.leftBegin()->getNodeID();
                SSNode left = getSSNode(leftID);
                // don't merge palindromic repeats / loops
                if (left.getFlag1())
                        break;
                if (left.numRightArcs() != 1)
                        break;
                nodeListq.push_front(leftID);
                left.setFlag1(true);
                curr = left;
        }

        if (nodeListq.size() <= 1)
                return result;  // only seed node

        // copy the deque into the vector
        copy(nodeListq.begin(), nodeListq.end(), std::back_inserter(result));
        return result;
}

void DBGraph::getConcatenations(std::vector<NodeChain>& reductions)
{
        reductions.clear();

        for (NodeID seedID = 1; seedID <= numNodes; seedID++) {
                NodeChain nc = getConcatenation(seedID);

                if (!nc.empty())
                        reductions.push_back(nc);
        }

        setAllFlags1(false);
}

void DBGraph::findLinearUniquePaths(NodeMap<int>& nodeMult,
                                    EdgeMap<int>& edgeMult)
{
        size_t numRemoved = 0;

        vector<NodeID> toConcatenate;

        for (NodeID sID = -numNodes; sID <= numNodes; sID++) {
                if (!nodeExists(sID))
                        continue;       // deleted node
                if (nodeMult[sID] != 1)
                        continue;       // not a unique node
                SSNode s = getSSNode(sID);
                if (s.numRightArcs() == 0)
                        continue;       // no right arcs

                vector<EdgeRep> toDelete;

                // check if there is a single right arc with mult = 1 pointing
                // to a node with multiplicity = 1
                int sumEdgeMult = 0; NodeID tID = 0;
                for (ArcIt it = s.rightBegin(); it != s.rightEnd(); it++) {
                        const EdgeRep er(sID, it->getNodeID());
                        sumEdgeMult += edgeMult[er];
                        if (edgeMult[er] == 1)
                                if (nodeMult[it->getNodeID()] == 1)
                                        tID = it->getNodeID();
                        if (edgeMult[er] == 0)
                                toDelete.push_back(er);
                }

                if ((sumEdgeMult != 1) || (tID == 0))
                        continue;

                SSNode t = getSSNode(tID);

                // check if there is a single left arc with mult = 1
                sumEdgeMult = 0;
                for (ArcIt it = t.leftBegin(); it != t.leftEnd(); it++) {
                        const EdgeRep er(it->getNodeID(), tID);
                        sumEdgeMult += edgeMult[er];
                        if (edgeMult[er] == 0)
                                toDelete.push_back(er);
                }

                if (sumEdgeMult != 1)
                        continue;

                for (const auto& er: toDelete)
                        removeEdges(toDelete, false);

                toConcatenate.push_back(sID);

                cout << "Concatenated around node " << sID << endl;

                numRemoved += toDelete.size();
        }

        for (NodeID id : toConcatenate) {
                if (!nodeExists(id))
                        continue;
                vector<NodeID> nodeList_v;
                concatenateAroundNode(id, nodeList_v);
        }

        //return numRemoved;
}
