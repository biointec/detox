/******************************************************************************
 *   Copyright (C) 2014 - 2026 Jan Fostier (jan.fostier@ugent.be)             *
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
#include "dbgraph.h"
#include "crfmult.h"
#include "alignment.h"
#include "correctgraph.h"
#include "pathinfo.h"
#include "nodetracker.h"

using namespace std;

// ============================================================================
// REMOVE NODES AND ARCS BASED ON COVERAGE
// ============================================================================

void DBGraph::removeCoverage(double covCutoff, size_t maxMargLength)
{
        size_t initNumValidArcs = numValidArcs;
        size_t initNumValidNodes = numValidNodes;

        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = getSSNode(id);

                if (!node.isValid())
                        continue;

                // check if there are arcs to delete
                vector<NodeID> arcToDelete;
                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++)
                        if (it->getCov() <= covCutoff)
                                arcToDelete.push_back(it->getNodeID());
                for (auto it : arcToDelete)
                        removeArc(id, it);

                arcToDelete.clear();
                for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++)
                        if (it->getCov() <= covCutoff)
                                arcToDelete.push_back(it->getNodeID());
                for (auto it : arcToDelete)
                        removeArc(it, id);

                if (node.getAvgCov() > covCutoff)
                        continue;

                if (node.getMarginalLength() > maxMargLength)
                        continue;

                removeNode(id);
        }

        cout << "Removing " << initNumValidNodes - numValidNodes << " nodes "
             << "and " << initNumValidArcs - numValidArcs << " edges "
             << "with coverage <= " << covCutoff << "\n";
        concatenateNodes();
        cout << "\tGraph has " << getNumValidNodes() << " nodes and "
             << getNumValidArcs() << " arcs" << endl;
}

// ============================================================================
// TIP CLIPPING
// ============================================================================

ParallelPath DBGraph::findSimilarParallelTip(NodeID tipID, NodeLength maxLen)
{
        // This routine finds, given tipID, the most similar path of the form:
        //       ->- orig1 ->- orig2 ->- tipID
        // srcID |
        //       ->- para1 ->- para2
        // where origX and tipID are nodes with a *single* incoming arc that
        // have collective length of at most maxLen.
        // The original and parallel path are disjoint.
        // If no path exist, or MAX_TIP_ITERATIONS is exceeded, returns {}

        ParallelPath ret;

        SSNode tip = getSSNode(tipID);          // tip node
        assert(tip.isValid());
        assert(tip.numLeftArcs() == 1);
        assert(tip.numRightArcs() == 0);

        if (tip.getMarginalLength() > maxLen)   // tip too long
                return ret;

        float bestScore = numeric_limits<float>::min();
        int counter = 0, minPathLen = 0;
        NWAligner aln(MAX_TIP_INDELS, 1, -1, -3);
        aln.reserveBanded(maxLen, 0);

        vector<NodeID> orig = { tipID };
        string X;

        while (counter <= MAX_TIP_ITERATIONS)
        {
                SSNode orig1 = getSSNode(orig.front());
                orig1.setFlag2(true);
                X = orig1.substr(Kmer::getK() - 1) + X;
                assert(!X.empty());

                // conduct a depth-first search
                vector<pair<NodeID, int>> stack;        // node-depth combo
                vector<NodeID> para;                    // parallel path
                vector<int> pathLen;                    // path length

                // add the initial nodes to the stack (right neigbors of src)
                NodeID srcID = orig1.leftBegin()->getNodeID();
                vector<pair<NodeID, double>> next = getRightNeigbors(srcID);
                for (const auto [rID, cov] : next)
                        stack.emplace_back(rID, 0);

                while (!stack.empty() && (++counter <= MAX_TIP_ITERATIONS))
                {
                        const auto [currID, depth] = stack.back();
                        stack.pop_back();
                        SSNode n = getSSNode(currID);

                        if (n.getFlag2())
                                continue;       // do not go through origX

                        // update parallel path and pathLen
                        para.resize(depth);
                        para.push_back(currID);
                        pathLen.resize(depth);

                        // produce alignment
                        int offsetY = pathLen.empty() ? 0 : pathLen.back();
                        size_t maxYLen = X.size() + MAX_TIP_INDELS - offsetY;
                        string Y = n.substr(Kmer::getK() - 1, maxYLen);

                        AlnRes2 ar = aln.alignBandedContd(X, Y, offsetY);
                        pathLen.push_back(ar.lenY);

                        float thisScore = (float)ar.score / (float)X.size();
                        float maxAtt = (float)ar.maxAtt / (float)X.size();

                        // candSol = at least the origX nodes are aligned AND
                        // cannot/should not extend alignment
                        // (dead end partial OR full-length alignment)
                        bool candSol = (ar.lenX >= minPathLen) &&
                                (n.numRightArcs() == 0) || (ar.lenX == X.size());

                        // betSol = better score of at least same length
                        // OR longer alignment at least the same score
                        bool betSol =
                                ((thisScore >  bestScore) && (ar.lenX >= ret.alnOrig)) ||
                                ((thisScore >= bestScore) && (ar.lenX >  ret.alnOrig));

                        if (candSol && betSol) {
                                ret = ParallelPath(orig, para, ar.lenX, ar.lenY);
                                bestScore = thisScore;
                        }

                        // impossible to improve: get out
                        if (maxAtt <= bestScore)
                                continue;

                        // add right neigbors to the stack, coverage high to low
                        next = getRightNeigbors(currID);        // low to high
                        for (const auto [rID, cov] : next)
                                stack.emplace_back(rID, depth + 1);
                }

                // try and go one node further upstream
                SSNode src = getSSNode(srcID);
                if (src.numLeftArcs() != 1)
                        break;          // multiple incoming arcs
                if ((X.size() + src.getMarginalLength()) > maxLen)
                        break;          // too long
                minPathLen += src.getMarginalLength();

                orig.insert(orig.begin(), srcID);
        }

        for (auto id: orig)
                getSSNode(id).setFlag2(false);  // reset flag2

        return ret;
}

size_t DBGraph::removeLowCovTips(double maxCov, size_t maxLen,
                                 const std::set<NodeRep>& toKeep,
                                 bool extendTips)
{
        // protect nodes from being deleted
        NodeMap<bool> flag;
        for (NodeRep nr : toKeep) {
                flag[nr] = true;
                getSSNode(nr).setFlag1(true);   // prevent concatenation
        }

        // generate a list of candidate tips
        vector<pair<double, NodeID>> lnp;

        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (!nodeExists(id))
                        continue;       // deleted node
                SSNode n = getSSNode(id);
                if ((n.numLeftArcs() != 1) || (n.numRightArcs() > 0))
                        continue;       // not a tip
                if (n.getMarginalLength() > maxLen)
                        continue;       // node too long
                if (n.getAvgCov() > maxCov)
                        continue;       // coverage too high
                lnp.emplace_back(n.getAvgCov(), id);
        }

        auto sortByCov = [] (pair<double, NodeID> const& s1,
                             pair<double, NodeID> const& s2) -> bool {
                if (s1.first < s2.first)
                        return true;
                return NodeRep(s1.second) < NodeRep(s2.second);
        };

        sort(lnp.begin(), lnp.end());   // sort coverage low to high

        size_t numRemoved = 0;
        for (const auto [ac, tipID] : lnp)
        {
                SSNode tip = getSSNode(tipID);
                if (!tip.isValid() || flag.find(tipID) != flag.end())
                        continue;       // tip is parallel to another tip

                // find parallel path with most similar sequence
                ParallelPath res = findSimilarParallelTip(tipID, maxLen);

                const vector<NodeID>& orig = res.orig;
                const vector<NodeID>& para = res.para;

                if (orig.empty() || para.empty())
                        continue;       // no parallel path found, not a tip

                NodeID srcID = getSSNode(orig.front()).leftBegin()->getNodeID();

                // reduce coverage of original path nodes
                NodeLength origLen = tip.getMarginalLength();
                for (size_t i = 0; i < orig.size() - 1; i++) {
                        SSNode o = getSSNode(orig[i]);
                        NodeLength ml = o.getMarginalLength();
                        origLen += ml;
                        o.setCov(max(0.0, o.getCov() - ac * ml));
                }

                // reduce coverage of original path arcs
                for (size_t i = 0; i < orig.size() - 1; i++) {
                        Arc& a = getArc(i == 0 ? srcID : orig[i-1], orig[i]);
                        Arc& r = getArc(-orig[i], i == 0 ? -srcID : -orig[i-1]);
                        Coverage newCov = max(0.0, a.getCov() - ac);
                        a.setCov(newCov);
                        r.setCov(newCov);
                }

                // check for the presence of an unaligned tail of the tip
                NodeLength tailLen = origLen - res.alnOrig;
                assert(tailLen <= tip.getMarginalLength());

                // WARNING: this might introduce existing k-mers. Check!
                if ((tailLen > 0) && extendTips) {
                        string tail = tip.substr(tip.length() - tailLen);
                        SSNode p = getSSNode(para.back());
                        assert(p.numRightArcs() == 0);
                        p.setSequence(p.getSequence() + tail);
                }

                // increase coverage of the nodes of the parallel path
                for (NodeID pID : para) {
                        SSNode p = getSSNode(pID);
                        NodeLength thisLen = min(p.getMarginalLength(), origLen);
                        origLen -= thisLen;
                        p.setCov(p.getCov() + ac * thisLen);
                        flag[pID] = true;
                }

                // transfer coverage the edges of the parallel path
                for (size_t i = 0; i < para.size(); i++)
                        incArcCov(i == 0 ? srcID : para[i-1], para[i], ac);

                // remove node and concatenate adjacent nodes
                removeNode(tipID, true);
                numRemoved++;
        }

        return numRemoved;
}

bool DBGraph::getTipsTree(NodeID srcID, set<NodeID>& tips, size_t maxLen)
{
        assert(getSSNode(srcID).isValid());   // sanity check

        tips.clear();

        vector<pair<NodeID, size_t>> stack;
        stack.emplace_back(srcID, 0);

        while (!stack.empty()) {
                auto [currID, currLen] = stack.back();
                stack.pop_back();

                SSNode n = getSSNode(currID);

                if (n.numLeftArcs() != 1)
                        return false;   // not a tree

                currLen += n.getMarginalLength();

                if (currLen > maxLen)
                        return false;   // path too long

                if (!tips.insert(currID).second)
                        return false;   // node already inserted (= loop)

                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++)
                        stack.emplace_back(it->getNodeID(), currLen);

        }

        return true;
}

size_t DBGraph::removeTipsCRF(size_t maxLen, CRFMult& myCRFMult,
                              CovModel& nodeModel, CovModel& edgeModel)
{
        size_t numRemoved = 0;

        for (NodeID sID = -numNodes; sID <= numNodes; sID++) {
                if (!nodeExists(sID))
                        continue;       // deleted node
                SSNode s = getSSNode(sID);
                if ((s.numRightArcs() != 2) || (s.getMarginalLength() < maxLen))
                        continue;       // not a big source node
                ArcIt it = s.rightBegin();
                NodeID tID = it->getNodeID();
                it++;
                NodeID pID = it->getNodeID();
                SSNode t = getSSNode(tID);
                SSNode p = getSSNode(pID);
                // identify who is the candidate tip and who the parallel path
                if (p.getAvgCov() < t.getAvgCov()) {
                        swap(tID, pID);
                        swap(t, p);
                }
                if (p.getMarginalLength() < maxLen)
                        continue;       // not a big parallel node
                set<NodeID> tips;
                if (!getTipsTree(tID, tips, maxLen))
                        continue;       // not a set of tips

                NodeMap<Multiplicity> nodeMult(tips.size() + 2);
                nodeMult[sID] = Multiplicity();
                nodeMult[pID] = Multiplicity();
                for (NodeID id : tips)
                        nodeMult[NodeRep(id)] = Multiplicity();

                EdgeMap<Multiplicity> edgeMult;
                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel, true);

                if (nodeMult[sID].getExpMult() == 0)
                        continue;
                if (nodeMult[pID].getExpMult() == 0)
                        continue;
                bool allTips = true;
                for (auto id : tips)
                        if (nodeMult[NodeRep(id)].getExpMult() != 0)
                                allTips = false;

                if (!allTips)
                        continue;

                for (NodeID id : tips) {
                        removeNode(id, true);
                        //cout << "Node " << id << " is a tip" << endl;
                        numRemoved++;
                }
        }

        /*for (NodeID tID = -numNodes; tID <= numNodes; tID++) {
                // tip node
                if (!nodeExists(tID))
                        continue;       // deleted node
                SSNode t = getSSNode(tID);
                if ((t.numLeftArcs() != 1) || (t.numRightArcs() > 0))
                        continue;       // not a tip
                if (t.getMarginalLength() > maxLen)
                        continue;       // node too long

                // source node
                NodeID sID = t.leftBegin()->getNodeID();
                SSNode s = getSSNode(sID);
                if (s.numRightArcs() != 2)
                        continue;       // not a canonical tip

                // parallel node
                NodeID pID = 0;
                for (ArcIt it = s.rightBegin(); it != s.rightEnd(); it++) {
                        if (it->getNodeID() == tID)
                                continue;
                        pID = it->getNodeID();
                }
                SSNode p = getSSNode(pID);

                // avoid self-loops, palindromic repeats
                if (NodeRep(pID) == NodeRep(sID))
                        continue;

                /*if ((s.getMarginalLength() < maxLen) &&
                    (p.getMarginalLength() < maxLen))
                        continue;

                // compute the CRF
                NodeMap<Multiplicity> nodeMult(3);
                nodeMult[tID] = Multiplicity();
                nodeMult[sID] = Multiplicity();
                nodeMult[pID] = Multiplicity();
                EdgeMap<Multiplicity> edgeMult;
                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel);

                if (nodeMult[tID].getExpMult() != 0)
                        continue;
                if (nodeMult[sID].getExpMult() != 1)
                        continue;
                if (nodeMult[pID].getExpMult() != 1)
                        continue;

                removeNode(tID, true);

                cout << "Node " << tID << " is a tip" << endl;
                numRemoved++;
        }*/

        return numRemoved;
}

size_t DBGraph::removeBubblesCRF(size_t maxLen, CRFMult& myCRFMult,
                                 CovModel& nodeModel, CovModel& edgeModel)
{
        size_t numRemoved = 0;

        for (NodeID bID = 1; bID <= numNodes; bID++) {
                // bubble node
                if (!nodeExists(bID))
                        continue;       // deleted node
                SSNode b = getSSNode(bID);
                if ((b.numLeftArcs() != 1) || (b.numRightArcs() != 1))
                        continue;       // not a bubble
                if (b.getMarginalLength() > maxLen)
                        continue;       // node too long

                // source node
                NodeID sID = b.leftBegin()->getNodeID();
                SSNode s = getSSNode(sID);
                if (s.numRightArcs() != 2)
                        continue;       // not a canonical bubble

                // parallel node
                NodeID pID = 0;
                for (ArcIt it = s.rightBegin(); it != s.rightEnd(); it++) {
                        if (it->getNodeID() == bID)
                                continue;
                        pID = it->getNodeID();
                }
                SSNode p = getSSNode(pID);
                if ((p.numLeftArcs() != 1) || (p.numRightArcs() != 1))
                        continue;       // not a canonical bubble

                // destination node
                NodeID dID = b.rightBegin()->getNodeID();
                if (dID != p.rightBegin()->getNodeID())
                        continue;       // not a canonical bubble
                SSNode d = getSSNode(dID);
                if (d.numLeftArcs() != 2)
                        continue;       // not a canonical bubble

                if ((s.getMarginalLength() < maxLen) &&
                    (d.getMarginalLength() < maxLen))
                        continue;

                // compute the CRF
                NodeMap<Multiplicity> nodeMult(4);
                nodeMult[bID] = Multiplicity();
                nodeMult[sID] = Multiplicity();
                nodeMult[pID] = Multiplicity();
                nodeMult[dID] = Multiplicity();
                EdgeMap<Multiplicity> edgeMult;
                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel, true);

                if (nodeMult[bID].getExpMult() != 0)
                        continue;
                if (nodeMult[sID].getExpMult() != 1)
                        continue;
                if (nodeMult[pID].getExpMult() != 1)
                        continue;
                if (nodeMult[dID].getExpMult() != 1)
                        continue;

                removeNode(bID, true);

                //cout << "Node " << bID << " is a bubble" << endl;
                numRemoved++;
        }

        return numRemoved;
}

size_t DBGraph::removeNodesCRF(size_t maxLen, CRFMult& myCRFMult,
                               CovModel& nodeModel, CovModel& edgeModel)
{
        size_t numRemoved = 0;

        for (NodeID sID = -numNodes; sID <= numNodes; sID++) {
                if (!nodeExists(sID))
                        continue;       // deleted node
                SSNode s = getSSNode(sID);
                if ((s.numRightArcs() != 2) || (s.getMarginalLength() < maxLen))
                        continue;       // not a big source node
                ArcIt it = s.rightBegin();
                NodeID tID = it->getNodeID();
                it++;
                NodeID pID = it->getNodeID();
                SSNode t = getSSNode(tID);
                SSNode p = getSSNode(pID);
                // identify who is the candidate tip and who the parallel path
                if (p.getAvgCov() < t.getAvgCov()) {
                        swap(tID, pID);
                        swap(t, p);
                }
                if (p.getMarginalLength() < maxLen)
                        continue;       // not a big parallel node
                if (t.getMarginalLength() > maxLen)
                        continue;       // not a big parallel node

                NodeMap<Multiplicity> nodeMult(3);
                nodeMult[sID] = Multiplicity();
                nodeMult[pID] = Multiplicity();
                nodeMult[tID] = Multiplicity();

                EdgeMap<Multiplicity> edgeMult;
                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel, true);

                if (nodeMult[sID].getExpMult() == 0)
                        continue;
                if (nodeMult[pID].getExpMult() == 0)
                        continue;
                if (nodeMult[tID].getExpMult() > 0)
                        continue;

                removeNode(tID, true);
                //cout << "Node " << id << " is a tip" << endl;
                numRemoved++;
        }

        return numRemoved;
}

ParallelPath DBGraph::findSimilarParallelPath(NodeID id)
{
        // This routine finds, given id, the most similar path of the form:
        //       ->- orig1 ->- orig2 ->- id ->- orig3 ->-
        // srcID |                                      | dstID
        //       ->- path1 ->- path2 ->-
        // where origX upstream of id are nodes with a *single* incoming and
        // where origX downstream of id are nodes with a *single* outgoing arc
        // collectively are smaller than MAX_BUBBLE_LEN.
        // The original and parallel path are completely disjoint.
        // If no path exist, or MAX_BUBBLE_ITERATIONS is exceeded, returns {}

        float bestScore = numeric_limits<float>::min();
        ParallelPath ret;

        SSNode noi = getSSNode(id);             // node-of-interest
        assert(noi.isValid());
        assert(noi.numLeftArcs() == 1);
        assert(noi.numRightArcs() == 1);

        if (noi.getMarginalLength() > MAX_BUBBLE_LEN)      // tip too long
                return ret;

        int counter = 0;
        NWAligner aln(MAX_BUBBLE_INDELS, 1, -1, -3);

        vector<NodeID> orig = { id };
        string X;

        // dstID = most downstream, non-branching node
        NodeID dstID = noi.rightBegin()->getNodeID();
        size_t rlen = 0;
        while (true) {
                SSNode dst = getSSNode(dstID);
                if (dst.numRightArcs() != 1)
                        break;
                if (rlen + dst.getMarginalLength() > MAX_BUBBLE_LEN)
                        break;
                orig.push_back(dstID);
                rlen += dst.getMarginalLength();
                X = X + dst.substr(Kmer::getK() - 1);
                dstID = dst.rightBegin()->getNodeID();
        }

        while (true)
        {
                SSNode orig1 = getSSNode(orig.front());
                X = orig1.getSequence().substr(Kmer::getK() - 1) + X;

                aln.reserveBanded(X.size(), 0);

                // conduct a depth-first search
                vector<pair<NodeID, int>> stack;        // node-depth combo
                vector<NodeID> path;                    // node path
                vector<int> pathLen;                    // path length

                // add the initial nodes to the stack (right neigbors of src)
                NodeID srcID = orig1.leftBegin()->getNodeID();
                vector<pair<NodeID, double>> next = getRightNeigbors(srcID);
                for (const auto [rID, cov] : next)
                        stack.emplace_back(rID, 0);

                while (!stack.empty() && (counter++ < MAX_BUBBLE_ITERATIONS)) {
                        const auto [currID, depth] = stack.back();
                        stack.pop_back();

                        if (NodeRep(currID) == NodeRep(orig.front()))
                                continue;       // do not go through orig1

                        // update path and pathLen
                        path.resize(depth);
                        path.push_back(currID);
                        pathLen.resize(depth);

                        // produce alignment
                        SSNode n = getSSNode(currID);
                        int offsetY = pathLen.empty() ? 0 : pathLen.back();
                        size_t maxYLen = X.size() + MAX_BUBBLE_INDELS - offsetY;
                        string Y = n.substr(Kmer::getK() - 1, maxYLen);
                        AlnRes2 ar = aln.alignBandedContd(X, Y, offsetY);
                        pathLen.push_back(ar.lenY);

                        float thisScore = (float)ar.score  / (float)X.size();
                        float maxAtt    = (float)ar.maxAtt / (float)X.size();

                        // candSol = destination reached
                        bool candSol = edgeExists(EdgeRep(currID, dstID));

                        // betSol = better score
                        bool betSol = thisScore > bestScore;

                        if (candSol && betSol) {
                                ret = ParallelPath(orig, path, ar.lenX, ar.lenY);
                                bestScore = thisScore;
                        }

                        // destination reached or impossible to improve: get out
                        if (edgeExists(EdgeRep(currID, dstID)) || (maxAtt <= bestScore))
                                continue;

                        // push right neigbors on the stack, coverage high to low
                        next = getRightNeigbors(currID);        // low to high
                        for (const auto [rID, cov] : next)
                                stack.emplace_back(rID, depth + 1);
                }

                // try and go one node further upstream
                SSNode src = getSSNode(srcID);
                if (src.numLeftArcs() != 1)
                        break;          // multiple incoming arcs
                if (src.leftBegin()->getNodeID() == srcID)
                        break;          // move to a *different* node
                if ((X.size() + src.getMarginalLength()) > MAX_BUBBLE_LEN)
                        break;          // too long

                orig.insert(orig.begin(), srcID);
        }

        return ret;
}

size_t DBGraph::removeLowCovBubbles(double maxCov, size_t maxLen,
                                    const set<NodeRep>& toKeep)
{
        // protect nodes from being deleted
        NodeMap<bool> flag;
        for (NodeRep nr : toKeep) {
                flag[nr] = true;
                getSSNode(nr).setFlag1(true);   // prevent concatentation
        }

        // generate a list of candidate bubbles, sorted according to coverage
        vector<pair<double, NodeID>> lnp;

        for (NodeID id = 1; id <= numNodes; id++) {
                if (!nodeExists(id))
                        continue;       // deleted node
                SSNode n = getSSNode(id);
                if ((n.numLeftArcs() != 1) || (n.numRightArcs() != 1))
                        continue;       // not a bubble
                if ((maxLen != 0) && (n.getMarginalLength() > maxLen))
                        continue;       // node too long
                if (n.getAvgCov() > maxCov)
                        continue;       // coverage too high
                lnp.emplace_back(n.getAvgCov(), id);
        }

        auto sortByCov = [] (pair<double, NodeID> const& s1,
                             pair<double, NodeID> const& s2) -> bool {
                if (s1.first < s2.first)
                        return true;
                return NodeRep(s1.second) < NodeRep(s2.second);
        };

        sort(lnp.begin(), lnp.end());   // sort coverage low to high

        size_t numRemoved = 0;
        for (const auto [ac, id] : lnp)
        {
                SSNode n = getSSNode(id);
                if (!n.isValid() || flag.find(id) != flag.end())
                        continue;       // node already serves as a parallel path

                // find parallel path with similar sequence and higher coverage
                ParallelPath res = findSimilarParallelPath(id);
                vector<NodeID>& orig = res.orig;
                vector<NodeID>& para = res.para;

                if (orig.empty())
                        continue;       // no parallel path found, not a tip

                // remove possible common tail
                while (!para.empty() && orig.back() == para.back()) {
                        orig.pop_back();
                        para.pop_back();
                        assert(!orig.empty());
                }

                NodeID srcID = getSSNode(orig.front()).leftBegin()->getNodeID();
                NodeID dstID = getSSNode(orig.back()).rightBegin()->getNodeID();

                // reduce coverage of original path nodes
                for (NodeID oID : orig) {
                        SSNode o = getSSNode(oID);
                        NodeLength ml = o.getMarginalLength();
                        o.setCov(max(0.0, o.getCov() - ac * ml));
                }

                // reduce coverage of original path arcs
                for (size_t i = 0; i <= orig.size(); i++) {
                        NodeID s = (i == 0) ? srcID : orig[i-1];
                        NodeID d = (i == orig.size()) ? dstID : orig[i];

                        Arc& a = getArc(s, d);
                        Arc& r = getArc(-d, -s);
                        Coverage newCov = max(0.0, a.getCov() - ac);
                        a.setCov(newCov);
                        r.setCov(newCov);
                }

                // increase coverage of the nodes of the parallel path
                for (NodeID pID : para) {
                        SSNode p = getSSNode(pID);
                        NodeLength ml = p.getMarginalLength();
                        p.setCov(p.getCov() + ac * ml);
                        flag[pID] = true;
                }

                // transfer coverage the edges of the parallel path
                for (size_t i = 0; i <= para.size(); i++) {
                        NodeID s = (i == 0) ? srcID : para[i-1];
                        NodeID d = (i == para.size()) ? dstID : para[i];

                        Arc& a = getArc(s, d);
                        Arc& r = getArc(-d, -s);
                        Coverage newCov = a.getCov() + ac;
                        a.setCov(newCov);
                        r.setCov(newCov);
                }

                // remove node and concatentate adjacent nodes
                removeNode(id, true);
                numRemoved++;
        }

        size_t numRemaining = 0;
        for (const auto [ac, id] : lnp)
        {
                SSNode n = getSSNode(id);
                if (!n.isValid() || n.getFlag1())
                        continue;       // tip is parallel to another tip
                numRemaining++;
        }

        //cout << "Number of nodes NOT removed: " << numRemaining << endl;

        setAllFlags1(false);    // reset all flags

        return numRemoved;
}

// ============================================================================
// BUBBLE CORRECTION
// ============================================================================

void DBGraph::removeNode(NodeID nodeID, bool concatenate)
{
        SSNode node = getSSNode(nodeID);
        if (!node.isValid())            // node already removed
                return;

        vector<NodeID> toConcatenate;
        for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++) {
                SSNode leftNode = getSSNode(it->getNodeID());
                if (abs(leftNode.getNodeID()) == abs(node.getNodeID()))
                        continue;       // skip self-loops or palindromic arcs
                toConcatenate.emplace_back(it->getNodeID());
                leftNode.deleteRightArc(node.getNodeID());
                numValidArcs--;
        }

        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                SSNode rightNode = getSSNode ( it->getNodeID() );
                if (abs(rightNode.getNodeID()) == abs(node.getNodeID()))
                        continue;       // skip self-loops or palindromic arcs
                toConcatenate.emplace_back(it->getNodeID());
                rightNode.deleteLeftArc(node.getNodeID());
                numValidArcs--;
        }

        numValidArcs -= node.numRightArcs();
        node.deleteAllRightArcs();
        numValidArcs -= node.numLeftArcs();
        node.deleteAllLeftArcs();

        numValidNodes--;
        node.invalidate();

        if (concatenate) {
                vector<NodeID> nodeList_v;
                for (auto id : toConcatenate)
                        concatenateAroundNode(id, nodeList_v);
        }
}


// ============================================================================
// CRF-BASED GRAPH CORRECTION
// ============================================================================

vector<NodeRep> DBGraph::getLowCovNodes(double threshold) const
{
        vector<NodeRep> nodes;
        nodes.reserve(getNumValidNodes());
        for (size_t i = 1; i <= numNodes; i++) {
                SSNode n = getSSNode(i);
                if (n.isValid() && (n.getAvgCov() <= threshold))
                        nodes.push_back(NodeRep(i));
        }

        return nodes;
}

void DBGraph::filterNodes(const vector<NodeRep>& nodeReps,
                          const vector<bool>& flowOK,
                          vector<NodeRep>& toRemove, vector<NodeRep>& toKeep,
                          int depth)
{
        // sort nodes according to ascending coverage
        vector<pair<double, NodeRep>> lnp;
        lnp.reserve(nodeReps.size());
        for (size_t i = 0; i < nodeReps.size(); i++) {
                if (!flowOK[i])
                        continue;
                NodeRep nr = nodeReps[i];
                Coverage cov = getSSNode(nr).getAvgCov();
                lnp.emplace_back(make_pair(cov, nr));
        }
        sort(lnp.begin(), lnp.end());

        // select nodes to remove
        NodeMap<bool> nodeFlag;
        for (const auto [ignore, nr] : lnp) {
                if (nodeFlag.find(nr) != nodeFlag.end()) {      // keep node
                        toKeep.push_back(nr);
                } else {                                        // remove node
                        toRemove.push_back(nr);
                        vector<NodeID> sgn; vector<EdgeID> sge;
                        getSubgraph(nr, sgn, sge, depth);
                        for (const auto& nbID : sgn)
                                nodeFlag[nbID] = true;
                }
        }
}

void DBGraph::filterEdges(const vector<EdgeRep>& edgeReps,
                          const vector<bool>& flowOK,
                          vector<EdgeRep>& toRemove, vector<EdgeRep>& toKeep,
                          int depth)
{
        // sort edges according to ascending coverage
        vector<pair<double, EdgeRep>> lnp;
        lnp.reserve(edgeReps.size());
        for (size_t i = 0; i < edgeReps.size(); i++) {
                if (!flowOK[i])
                        continue;
                EdgeRep er = edgeReps[i];
                Coverage cov = getArc(er).getCov();
                lnp.emplace_back(make_pair(cov, er));
        }
        sort(lnp.begin(), lnp.end());

        // select edges to remove
        EdgeMap<bool> edgeFlag;
        for (const auto [ignore, er] : lnp) {
                if (edgeFlag.find(er) != edgeFlag.end()) {      // keep edge
                        toKeep.push_back(er);
                } else {                                        // remove edge
                        toRemove.push_back(er);
                        vector<NodeID> sgn; vector<EdgeID> sge;
                        getSubgraph(er.getSrcID(), sgn, sge, depth);
                        for (const auto& nbID : sge)
                                edgeFlag[nbID] = true;
                }
        }
}

void DBGraph::removeLowCovNodesFlow(double covCutoff, CRFMult& myCRFMult,
                                    CovModel& nodeModel, CovModel& edgeModel)
{
        vector<NodeRep> nodeReps;
        nodeReps.reserve(getNumValidNodes());
        for (size_t i = 1; i <= numNodes; i++) {
                SSNode n = getSSNode(i);
                if (n.isValid() && (n.getAvgCov() <= covCutoff) &&
                   (nodeModel.getExpMult(n.getAvgCov()) == 0))
                        nodeReps.emplace_back(NodeRep(i));
        }
        cout << "Selected " << nodeReps.size()
             << " nodes with coverage <= " << covCutoff << endl;

        while (!nodeReps.empty()) {
                vector<bool> flowOK;
                myCRFMult.checkFlow(nodeReps, flowOK, nodeModel, edgeModel);

                vector<NodeRep> toRemove, toKeep;
                filterNodes(nodeReps, flowOK, toRemove, toKeep,
                            myCRFMult.getMaxGraphDepth());

                cout << "\tRemoving " << toRemove.size() << " nodes" << endl;
                removeNodes(toRemove, true);
                cout << "\tGraph has " << getNumValidNodes()
                     << " nodes and " << getNumValidArcs() << " arcs\n";

                // select nodes to reconsider
                nodeReps.clear();
                for (NodeID id : toKeep) {
                        SSNode n = getSSNode(id);
                        if (n.isValid() && (n.getAvgCov() <= covCutoff))
                                nodeReps.push_back(id);
                }
        }
}

void DBGraph::removeLowCovEdgesFlow(double covCutoff, CRFMult& myCRFMult,
                                    CovModel& nodeModel, CovModel& edgeModel)
{
        vector<EdgeRep> edgeReps;
        edgeReps.reserve(getNumValidArcs() / 2);   // roughly half of the nodes
        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;
                SSNode n = getSSNode(id);
                if (!n.isValid())
                        continue;
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        if (id > -it->getNodeID() || it->getCov() > covCutoff)
                                continue;
                        if (edgeModel.getExpMult(it->getCov()) != 0)
                                continue;
                        edgeReps.push_back(EdgeRep(id, it->getNodeID()));
                }
        }

        cout << "Selected " << edgeReps.size()
             << " edges with coverage <= " << covCutoff << endl;

        while (!edgeReps.empty()) {
                vector<bool> flowOK;
                myCRFMult.checkFlow(edgeReps, flowOK, nodeModel, edgeModel);

                vector<EdgeRep> toRemove, toKeep;
                filterEdges(edgeReps, flowOK, toRemove, toKeep,
                            myCRFMult.getMaxGraphDepth());

                cout << "\tRemoving " << toRemove.size() << " edges" << endl;
                removeEdges(toRemove, true);
                cout << "\tGraph has " << getNumValidNodes()
                     << " nodes and " << getNumValidArcs() << " arcs\n";

                // select edges to reconsider
                edgeReps.clear();
                for (EdgeRep er : toKeep) {
                        if (edgeExists(er) && (getArc(er).getCov() <= covCutoff))
                                edgeReps.push_back(er);
                }
        }
}

void DBGraph::filterNodes(const NodeMap<Multiplicity>& nodeMult,
                          vector<NodeRep>& toRemove, vector<NodeRep>& toKeep,
                          int depth)
{
        // sort nodes according to descending log odd-ratio
        vector<pair<double, NodeRep>> lnp;
        lnp.reserve(nodeMult.size());
        for (const auto& it : nodeMult) {
                if (it.second.getExpMult() != 0)
                        continue;
                const double& LOR = it.second.getExpMultLogOR();
                lnp.emplace_back(make_pair(LOR, it.first));
        }
        sort(lnp.begin(), lnp.end());
        reverse(lnp.begin(), lnp.end());

        // select nodes to remove
        for (const auto& e : lnp) {
                NodeRep id = e.second;
                if (getSSNode(id).getFlag1()) {         // keep node
                        toKeep.push_back(id);
                } else {                                // remove node
                        toRemove.push_back(id);
                        vector<NodeID> sgn; vector<EdgeID> sge;
                        getSubgraph(id, sgn, sge, depth);
                        for (const auto& nbID : sgn)
                                getSSNode(nbID).setFlag1(true);
                }
        }

        setAllFlags1(false);
}

void DBGraph::removeLowCovNodesCRF(double covCutoff, CRFMult& myCRFMult,
                                   CovModel& nodeModel, CovModel& edgeModel)
{
        vector<NodeRep> nodeReps = getLowCovNodes(covCutoff);
        cout << "Selected " << nodeReps.size()
             << " nodes with coverage <= " << covCutoff << endl;
        while (!nodeReps.empty()) {
                NodeMap<Multiplicity> nodeMult(nodeReps.size());
                for (size_t i = 0; i < nodeReps.size(); i++)
                        nodeMult[nodeReps[i]] = Multiplicity();
                EdgeMap<Multiplicity> edgeMult;
                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel);

                vector<NodeRep> toRemove, toKeep;
                filterNodes(nodeMult, toRemove, toKeep,
                            myCRFMult.getMaxGraphDepth());

                cout << "\tRemoving " << toRemove.size() << " nodes" << endl;
                removeNodes(toRemove, true);
                cout << "\tGraph has " << getNumValidNodes()
                     << " nodes and " << getNumValidArcs() << " arcs\n";

                // select nodes to reconsider
                nodeReps.clear();
                for (NodeID id : toKeep) {
                        SSNode n = getSSNode(id);
                        if (n.isValid() && (n.getAvgCov() <= covCutoff))
                                nodeReps.push_back(id);
                }
        }
}

void DBGraph::filterEdges(const EdgeMap<Multiplicity>& edgeMult,
                          vector<EdgeRep>& toRemove, vector<EdgeRep>& toKeep,
                          int depth)
{
        // sort edges according to descending log odd-ratio
        vector<pair<double, EdgeRep>> lnp;
        for (const auto& it : edgeMult) {
                if (it.second.getExpMult() != 0)
                        continue;
                const double& LOR = it.second.getExpMultLogOR();
                lnp.push_back(make_pair(LOR, it.first));
        }
        sort(lnp.begin(), lnp.end());
        reverse(lnp.begin(), lnp.end());

        // select edges to remove
        EdgeMap<bool> edgeFlag;
        for (const auto& e : lnp) {
                EdgeRep id = e.second;
                if (edgeFlag.find(id) != edgeFlag.end()) {      // keep edge
                        toKeep.push_back(id);
                } else {                                        // remove edge
                        toRemove.push_back(id);
                        vector<NodeID> sgn; vector<EdgeID> sge;
                        getSubgraph(id.getSrcID(), sgn, sge, depth);
                        for (const auto& nbID : sge)
                                edgeFlag[nbID] = true;
                }
        }
}

void DBGraph::removeLowCovEdgesCRF(double covCutoff, CRFMult& myCRFMult,
                                   CovModel& nodeModel, CovModel& edgeModel)
{
        vector<EdgeRep> edgeReps = getLowCovEdges(covCutoff);
        cout << "Selected " << edgeReps.size()
             << " edges with coverage <= " << covCutoff << endl;
        while (!edgeReps.empty()) {
                EdgeMap<Multiplicity> edgeMult(edgeReps.size());
                for (size_t i = 0; i < edgeReps.size(); i++)
                        edgeMult[edgeReps[i]] = Multiplicity();
                NodeMap<Multiplicity> nodeMult;
                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel);

                vector<EdgeRep> toRemove, toKeep;
                filterEdges(edgeMult, toRemove, toKeep,
                            myCRFMult.getMaxGraphDepth());

                cout << "\tRemoving " << toRemove.size() << " edges" << endl;
                removeEdges(toRemove, true);
                cout << "\tGraph has " << getNumValidNodes()
                     << " nodes and " << getNumValidArcs() << " arcs\n";

                // select edges to reconsider
                edgeReps.clear();
                for (EdgeRep er : toKeep) {
                        if (edgeExists(er) && (getArc(er).getCov() <= covCutoff))
                                edgeReps.push_back(er);
                }
        }
}

// ============================================================================
// GENERIC ROUTINES FOR GRAPH CORRECTION
// ============================================================================

void DBGraph::removeArc(NodeID leftID, NodeID rightID, bool concatenate)
{
        if (getSSNode(leftID).deleteRightArc(rightID))
                numValidArcs--;

        if (getSSNode(rightID).deleteLeftArc(leftID))
                numValidArcs--;

        if (concatenate) {
                vector<NodeID> v;
                concatenateAroundNode(leftID, v);
                concatenateAroundNode(rightID, v);
        }
}

void DBGraph::removeRightArcs(NodeID nodeID)
{
        SSNode n = getSSNode(nodeID);

        vector<NodeID> rNeighbors;
        for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++)
                rNeighbors.push_back(it->getNodeID());

        for (NodeID rID : rNeighbors)
                removeArc(nodeID, rID);
}

void DBGraph::removeLeftArcs(NodeID nodeID)
{
        SSNode n = getSSNode(nodeID);

        vector<NodeID> lNeighbors;
        for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++)
                lNeighbors.push_back(it->getNodeID());

        for (NodeID lID : lNeighbors)
                removeArc(lID, nodeID);
}

void DBGraph::removeCoverageNodes(double covCutoff, size_t maxMargLength)
{
        size_t initNumValidArcs = numValidArcs;
        size_t initNumValidNodes = numValidNodes;

        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = getSSNode(id);

                if (!node.isValid())
                        continue;

                // check if there are arcs to delete
                vector<NodeID> arcToDelete;
                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++)
                        if (it->getCov() <= 0.0)
                                arcToDelete.push_back(it->getNodeID());
                for (auto it : arcToDelete)
                        removeArc(id, it);

                arcToDelete.clear();
                for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++)
                        if (it->getCov() <= 0.0)
                                arcToDelete.push_back(it->getNodeID());
                for (auto it : arcToDelete)
                        removeArc(it, id);

                if (node.getAvgCov() > covCutoff)
                        continue;

                if (node.getMarginalLength() > maxMargLength)
                        continue;

                removeNode(id);
        }

        cout << "\tRemoved " << initNumValidNodes - numValidNodes << " nodes\n";
        cout << "\tRemoved " << initNumValidArcs - numValidArcs << " arcs\n";
}

string DBGraph::convertNodesToString(const NodeChain& nc) const
{
        string ret;

        if (nc.empty())
                return ret;

        // handle circular node sequences
        size_t numNodes = ((nc.size() > 1) && (nc.front() == nc.back())) ?
                nc.size() - 1 : nc.size();

        size_t size = 0;
        for (size_t i = 0; i < numNodes; i++)
                size += getSSNode(nc[i]).length();
        size -= (numNodes - 1) * (Kmer::getK() - 1);

        ret = getSSNode(nc[0]).getSequence();
        for (size_t i = 1; i < numNodes; i++)
                ret.append(getSSNode(nc[i]).getSequence().substr(Kmer::getK()-1));

        assert(size == ret.size());
        return ret;
}

string DBGraph::convertNodesToShortestString(const NodeChain& nc) const
{
        string ret;

        if (nc.empty())
                return ret;     // return empty string

        ret = getSSNode(nc.front()).getRightKmer().str();
        if (nc.size() == 1)
                return ret;     // only a single node

        for (size_t i = 1; i < nc.size() - 1; i++)
                ret.append(getSSNode(nc[i]).getSequence().substr(Kmer::getK()-1));
        ret.append(getSSNode(nc.back()).getLeftKmer().str().substr(Kmer::getK()-1));

        return ret;
}

NodeChain DBGraph::convertShortestStringToNodes(NodeID srcID,
                                                const string& str) const
{
        NodeChain ret;

        if (str.size() < Kmer::getK())
                return ret;     // return empty node chain

        if (str.substr(0, Kmer::getK()) != getSSNode(srcID).getRightKmer().str())
                return ret;     // first k-mer != last k-mer of srcID

        ret.emplace_back(srcID);
        for (size_t i = Kmer::getK(); i < str.size(); ) {
                NodeID r = getSSNode(ret.back()).getRightArcNodeID(str[i]);
                if (r == 0)
                        return ret;
                ret.emplace_back(r);
                i += getSSNode(r).getMarginalLength();
        }

        return ret;
}

ConsNodeChain DBGraph::convertCNC(const ConsNodeChain& thisCNC,
                                  const DBGraph& newGraph,
                                  const BiNodeMap<NodeID>& thisToNew) const
{
        ConsNodeChain newCNC;

        NodeLength newLen = 0, newOffset = 0;
        for (auto [thisNodeID, thisWeight] : thisCNC) {
                NodeLength thisOffset = 0;
                NodeLength thisLen = getSSNode(thisNodeID).getMarginalLength();

                if (newOffset < newLen) {       // newNode is longer than origNode
                        NodeID newNodeID = newCNC.back().first;
                        char c = newGraph.getSSNode(newNodeID).getNucleotide(newOffset + Kmer::getK() - 1);
                        if (getSSNode(thisNodeID).peekNucleotideMarginalLeft() != c)
                                return newCNC;
                        newCNC.back().second = min<int>(newCNC.back().second, thisWeight);
                        thisOffset = min<NodeLength>(thisLen, newLen - newOffset);
                        newOffset += thisLen;
                }

                while (thisOffset < thisLen) {  // fill origNode with newNodes
                        char c = getSSNode(thisNodeID).getNucleotide(thisOffset + Kmer::getK() - 1);
                        NodeID newNodeID = (newCNC.empty()) ? get(thisToNew, thisNodeID) :
                                newGraph.getSSNode(newCNC.back().first).getRightArcNodeID(c);
                        if (newNodeID == 0)
                                return newCNC;
                        newLen = newGraph.getSSNode(newNodeID).getMarginalLength();
                        newCNC.emplace_back(newNodeID, thisWeight);
                        newOffset = min<NodeLength>(newLen, thisLen - thisOffset);
                        thisOffset += newLen;
                }
        }

        return newCNC;
}

NodeChain DBGraph::convertNC(const NodeChain& thisNC,
                             const DBGraph& newGraph,
                             const BiNodeMap<NodePosPair>& origToSub) const
{
        NodeChain newNC;

        NodeLength newLen = 0, newOffset = 0;
        for (size_t i = 0; i < thisNC.size(); i++) {
                NodeID thisNodeID = thisNC[i];
                SSNode thisNode = getSSNode(thisNodeID);
                NodeLength thisOffset = 0;
                NodeLength thisLen = thisNode.getMarginalLength();

                if (newOffset < newLen) {       // newNode is longer than origNode
                        NodeID newNodeID = newNC.back();
                        char c = newGraph.getSSNode(newNodeID).getNucleotide(newOffset + Kmer::getK() - 1);
                        if (thisNode.peekNucleotideMarginalLeft() != c)
                                return newNC;
                        thisOffset = min<NodeLength>(thisLen, newLen - newOffset);
                        newOffset += thisLen;
                }

                while (thisOffset < thisLen) {  // fill origNode with newNode(s)
                        NodeID newNodeID;
                        if (i == 0) {
                                NodePosPair npp = get(origToSub, thisNodeID);
                                newNodeID = npp.first;
                                newOffset = npp.second;
                        } else {
                                char c = thisNode.getNucleotide(thisOffset + Kmer::getK() - 1);
                                newNodeID = newGraph.getSSNode(newNC.back()).getRightArcNodeID(c);
                                newOffset = 0;
                        }

                        if (newNodeID == 0)
                                return newNC;
                        newLen = newGraph.getSSNode(newNodeID).getMarginalLength();
                        newNC.emplace_back(newNodeID);
                        newOffset = min<NodeLength>(newLen, thisLen - thisOffset + newOffset);
                        thisOffset += newLen;
                }
        }

        return newNC;
}

void DBGraph::concatenateAroundNode(NodeID seedID, vector<NodeID>& nodeListv)
{
        nodeListv.clear();

        SSNode seed = getSSNode(seedID);
        if (!seed.isValid())
                return;

        deque<NodeID> nodeListq;
        nodeListq.push_back(seedID);

        if (seed.getFlag1())
                return;         // flag 1 protects a node from concatenation

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

        // copy the deque into the vector
        copy(nodeListq.begin(), nodeListq.end(), std::back_inserter(nodeListv));

        // reset the flags to false
        for (const auto& it : nodeListq)
                getSSNode(it).setFlag1(false);

        // if no linear path was found, continue
        if (nodeListq.size() == 1)
                return;

        // concatenate the path
        NodeID frontID = nodeListq.front();
        SSNode front = getSSNode(frontID);
        NodeID backID = nodeListq.back();
        SSNode back = getSSNode(backID);

        front.deleteAllRightArcs();
        front.inheritRightArcs(back);

        string str = convertNodesToString(nodeListv);

        Coverage newCov = front.getCov();
        for (size_t i = 1; i < nodeListq.size(); i++) {
                SSNode n = getSSNode(nodeListq[i]);
                newCov += n.getCov();

                removeNode(nodeListq[i]);
        }

        front.setCov(newCov);
        front.setSequence(str);
}

bool DBGraph::concatenateNodes(const set<NodeRep>& toKeep)
{
        // protect nodes from being deleted
        for (NodeRep nr : toKeep)
                getSSNode(nr).setFlag1(true);

        size_t numConcatenations = 0;

        for (NodeID seedID = 1; seedID <= numNodes; seedID++) {
                vector<NodeID> concatenation;
                concatenateAroundNode(seedID, concatenation);

                if (!concatenation.empty())
                        numConcatenations += concatenation.size() - 1;
        }

        // undo flags
        for (NodeRep nr : toKeep)
                getSSNode(nr).setFlag1(false);

        //cout << "Concatenated " << numConcatenations << endl;

        return (numConcatenations > 0);
}

void DBGraph::removeNodes(const std::vector<NodeRep>& nodes, bool concatenate)
{
        for (const auto& el : nodes)
                removeNode(el, concatenate);
}

void DBGraph::removeEdges(const std::vector<EdgeRep>& edges, bool concatenate)
{
        for (const auto& el : edges)
                removeArc(el.getSrcID(), el.getDstID(), concatenate);
}

void DBGraph::getPERReductions(const LibraryContainer& library,
                            NodeMap<Multiplicity>& nodeMap,
                            EdgeMap<Multiplicity>& edgeMap,
                            vector<NodeChain>& red)
{
        vector<NodeChain> readsL, readsR;
        size_t j = 0;
        for (const Library& lib : library) {
                const auto [fn1, fn2] = lib.getAlnFilename();

                ifstream ifs1(fn1);
                ifstream ifs2(fn2);
                string line1, line2;
                while (getline(ifs1, line1) && getline(ifs2, line2)) {
                        NodeChain L, R;

                        istringstream iss(line1);
                        NodeID nodeID;
                        while (iss >> nodeID)
                                L.push_back(nodeID);

                        iss = istringstream(line2);
                        while (iss >> nodeID)
                                R.push_back(nodeID);

                        // FIXME: handle other orientations?
                        R = R.getRevCompl();

                        // we want proper paired-end reads
                        if (L.empty() || R.empty())
                                continue;

                        readsL.push_back(L);
                        readsR.push_back(R);
                }
        }

        NodeID noi = 2;

        map<pair<NodeID, NodeID>, int> uniquePairs;
        for (size_t i = 0; i < readsL.size(); i++) {
                NodeChain& L = readsL[i];
                NodeChain& R = readsR[i];

                // both L and R do not contain the noi
                if (!L.contains(NodeRep(noi)) && !R.contains(NodeRep(noi)))
                        continue;

                // if both L and R contain the noi, L and R should be mergeable
                // get out
                if (L.contains(NodeRep(noi)) && R.contains(NodeRep(noi)))
                        continue;

                // L is upstream from noi, get out (FIXME: information from R ??)
                if (R.contains(noi))
                        continue;
                if (L.contains(-noi))
                        continue;

                if (R.contains(-noi)) {
                        swap(L, R);
                        L = L.getRevCompl();
                        R = R.getRevCompl();
                }

                cout << L << "\t" << R << endl;

                // get the rightmost unique node of L
                NodeID nodeL = 0;
                for (size_t i = L.size(); i--> 0; )
                        if (nodeMap[L[i]].getExpMult() == 1) {
                                nodeL = L[i];
                                break;
                        }

                // get the leftmost unique node of R
                NodeID nodeR = 0;
                for (size_t i = 0; i < R.size(); i++)
                        if (nodeMap[R[i]].getExpMult() == 1) {
                                nodeR = R[i];
                                break;
                        }

                if (nodeL == 0 || nodeR == 0)
                        continue;

                uniquePairs[make_pair(nodeL, nodeR)]++;
        }

        for (auto it : uniquePairs) {
                cout << it.first.first << "," << it.first.second << "\t" << it.second << endl;
        }
}

void DBGraph::adjustMult(const vector<NodeChain>& reductions,
                         const vector<NodeChain>& looseEnds,
                         NodeMap<int>& estNodeMult, EdgeMap<int>& estEdgeMult,
                         NodeMap<int>& redNodeMult, EdgeMap<int>& redEdgeMult)
{
        redNodeMult.clear();
        redEdgeMult.clear();

        // count the multiplicities of the nodes/edges of reductions
        for (const auto& r : reductions) {
                for (size_t i = 0; i < r.size(); i++)
                        redNodeMult[NodeRep(r[i])]++;
                // correct for circular reduction (don't count twice)
                if ((r.size() > 1) && (r.front() == r.back()))
                        redNodeMult[NodeRep(r.front())]--;
                for (size_t i = 1; i < r.size(); i++)
                        redEdgeMult[EdgeRep(r[i-1], r[i])]++;
        }

        // make sure the first/last node of a reduction have multiplicity = 1
        for (const auto& r : reductions) {
                assert(redNodeMult[r.front()] == 1);
                assert(redNodeMult[r.back()] == 1);
        }

        // take into account the loose ends (they might be resolved further)
        NodeMap<int> leNodeMult; EdgeMap<int> leEdgeMult;
        for (const auto& r : looseEnds) {
                NodeMap<int> thisNodeMult; EdgeMap<int> thisEdgeMult;
                // nodes (exclude start node !!!)
                for (size_t i = 1; i < r.size(); i++)
                        thisNodeMult[NodeRep(r[i])]++;
                for (const auto& [nr, m] : thisNodeMult)
                        leNodeMult[nr] = max(m, leNodeMult[nr]);
                // edges
                for (size_t i = 1; i < r.size(); i++)
                        thisEdgeMult[EdgeRep(r[i-1], r[i])]++;
                for (const auto& [er, m] : thisEdgeMult)
                        leEdgeMult[er] = max(m, leEdgeMult[er]);
        }

        // if the estimated multiplicity is lower than the number of times
        // it is effectively used, raise it accordingly
        for (const auto& [nr, m] : redNodeMult)
                estNodeMult[nr] = max(estNodeMult[nr], m + get(leNodeMult, nr));
        for (auto& [er, m] : redEdgeMult)
                estEdgeMult[er] = max(estEdgeMult[er], m + get(leEdgeMult, er));

        // make sure the first/last node of a reduction have multiplicity = 1
        for (const auto& r : reductions) {
                estNodeMult[r.front()] = 1;
                estNodeMult[r.back()] = 1;
        }
}

void DBGraph::applyReductions(const std::vector<NodeChain>& reductions,
                              NodeMap<int>& nodeMult, EdgeMap<int>& edgeMult)
{
        for (const auto& r : reductions) {
                assert(r.size() > 1);
                SSNode nodeFront = getSSNode(r.front());
                SSNode nodeBack = getSSNode(r.back());

                // adjust the coverage and multiplicities of the internal nodes
                Coverage newCov = (r.front() == r.back()) ? nodeFront.getCov() :
                        nodeFront.getCov() + nodeBack.getCov();
                for (size_t i = 1; i < r.size() - 1; i++) {
                        SSNode n = getSSNode(r[i]);
                        Coverage oldCov = n.getCov();
                        int oldMult = nodeMult[r[i]];
                        assert(oldMult > 0);

                        // divide the coverage among new and old node
                        float f = float(oldMult - 1) / float(oldMult);
                        n.setCov(f * oldCov);
                        newCov += (1.0-f) * oldCov;

                        nodeMult[r[i]]--;
                }
                nodeFront.setCov(newCov);

                // adjust the coverage and multiplicities of the internal edges
                for (size_t i = 1; i < r.size() - 2; i++) {
                        EdgeRep e(r[i], r[i+1]);
                        Coverage oldCov = getArc(e).getCov();
                        int oldMult = edgeMult[e];
                        assert(oldMult > 0);

                        // reduce the coverage proportionally
                        float f = float(oldMult - 1) / float(oldMult);
                        getArc(r[i], r[i+1]).setCov(f * oldCov);
                        getArc(-r[i+1], -r[i]).setCov(f * oldCov);

                        edgeMult[e]--;
                }

                // set the concatenated sequence to node A
                string str = convertNodesToString(r);
                nodeFront.setSequence(str);

                // remove arcs
                removeRightArcs(r.front());
                removeLeftArcs(r.back());

                // transfer nodeC's right arcs to nodeA
                nodeFront.inheritRightArcs(nodeBack);

                // remove nodeC that now exists as nodeA
                if (r.front() != r.back())
                        removeNode(r.back());

                // remove edges with zero multiplicity
                for (size_t i = 1; i < r.size() - 2; i++)
                        if (edgeMult[EdgeRep(r[i], r[i+1])] == 0)
                                removeArc(r[i], r[i+1]);

                // remove nodes with zero multiplicity
                for (size_t i = 1; i < r.size() - 1; i++) {
                        if (nodeMult[r[i]] == 0)
                                removeNode(r[i]);
                }
        }
}

void DBGraph::glueTips(NodeID leftID, NodeID rightID)
{
        SSNode left = getSSNode(leftID);
        SSNode right = getSSNode(rightID);

        if (!left.isValid() || !right.isValid()) {
                cout << "Warning: skipped glueing of deleted nodes..." << endl;
                return;
        }

        string s1 = left.getSequence().substr(left.length() - MIN_TIP_MERGE_LEN);
        string s2 = right.getSequence().substr(0, MIN_TIP_MERGE_LEN);

        NWAligner aligner(3, 1, -1, -3);
        AlnRes res = aligner.overlapAln(s1, s2);

        string str = left.getSequence() + right.getSequence().substr(res.s2len);

        left.setSequence(str);
        left.inheritRightArcs(right);
        left.setCov(left.getCov() + right.getCov());

        removeNode(rightID);
}

void DBGraph::glueTips(LibraryContainer& libraries)
{
        BiNodeMap<string> tips;         // tips with their terminal string
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode n = getSSNode(id);
                if (!n.isValid())
                        continue;
                if (n.getMarginalLength() < MIN_TIP_MERGE_LEN)
                        continue;       // too short to merge
                if (n.numRightArcs() == 0)
                        tips[id] = n.getSequence().substr(n.length() - MIN_TIP_MERGE_LEN);
                n = getSSNode(-id);
                if (n.numRightArcs() == 0)
                        tips[-id] = n.getSequence().substr(n.length() - MIN_TIP_MERGE_LEN);
        }

        BiNodeMap<bool> noi;             // indicate the unique nodes
        for (const auto& [id, ignore] : tips)
                noi[id] = true;

        // build the path info from the aligned reads
        BiNodeMap<PathInfo> pathInfo;
        addPERInfo(libraries, noi, pathInfo);

        // make all pairwise overlaps and store, per tip, the best matching tip
        NWAligner aligner(3, 1, -1, -3);
        map<NodeID, pair<NodeID, int> > pairedNodes;
        for (const auto& [id1, str1] : tips) {
                BiNodeMap<int> perDst = pathInfo[id1].getUniquePerDst(noi);

                // find the best destination tip
                NodeID id2 = 0; int max = 0, secondMax = 0;
                for (const auto [dst, count] : perDst) {
                        if (NodeRep(id1) == NodeRep(dst))
                                continue;       // do not connect to yourself
                        if (count > max) {
                                id2 = dst;
                                secondMax = max;
                                max = count;
                        } else if (count > secondMax) {
                                secondMax = count;
                        }
                }

                if (id2 == 0)
                        continue;       // no target found

                if (secondMax * 4 > max)
                        continue;       // maximim conflict degree exceeded

                if (max < 2)
                        continue;       // FIXME: parameters

                if (tips.find(-id2) == tips.end())
                        continue;
                string str2 = Nucleotide::getRevCompl(tips[-id2]);
                AlnRes res = aligner.overlapAln(str1, str2);

                //cout << res.score << endl;

                // skip poor alignments
                if (res.score < int(Kmer::getK() / 2))
                        continue;

                if (res.score > pairedNodes[id1].second)
                        pairedNodes[id1] = make_pair(id2, res.score);

                //cout << id1 << " - " << -id2 << ": " << res.score << endl;
        }

        // find mutually best pairs of tips and glue them together
        int counter = 0;
        for (const auto& it : pairedNodes) {
                NodeID leftID = it.first;
                NodeID rightID = it.second.first;

                // stay out of trouble
                if (abs(leftID) == abs(rightID))
                        continue;
                // only keep one representation
                if (leftID > -rightID)
                        continue;
                // check if mutually best
                if (pairedNodes[-rightID].first != -leftID)
                        continue;

                //cout << "Glueing " << leftID << " to " << rightID << endl;
                /*{
                SSNode n = getSSNode(leftID);
                cout << "Node: " << n.getNodeID() << ", valid = "
                             << (n.isValid() ? "true " : "false ") << ", marg. len = "
                             <<  n.getMarginalLength() << ", avg. cov. = " << n.getAvgCov()
                             << "; arcs (l/r): " << (int)n.numLeftArcs() << "/" << (int)n.numRightArcs()
                             << "; mult (true/est): " << endl;
                             // << trueNodeMult[n.getNodeID()] << "/"
                             // << nodeMult[n.getNodeID()].getExpMult() << endl;
                }
                {
                SSNode n = getSSNode(rightID);
                cout << "Node: " << n.getNodeID() << ", valid = "
                             << (n.isValid() ? "true " : "false ") << ", marg. len = "
                             <<  n.getMarginalLength() << ", avg. cov. = " << n.getAvgCov()
                             << "; arcs (l/r): " << (int)n.numLeftArcs() << "/" << (int)n.numRightArcs()
                             << "; mult (true/est): " << endl;
                             // << trueNodeMult[n.getNodeID()] << "/"
                             // << nodeMult[n.getNodeID()].getExpMult() << endl;
                }*/

                glueTips(leftID, rightID);
                counter++;
        }

        cout << "Connected " << counter << " pairs of tips" << endl;
}

void DBGraph::getSubgraph(NodeID srcID, int maxLen, int maxNodes,
                          map<NodeID, int>& nodeDist)
{
        nodeDist.clear();

        priority_queue<NodeDepth, vector<NodeDepth>, NodeDepthComp> pq;
        pq.push(NodeDepth(srcID, 0));

        while (!pq.empty()) {
                // get and erase the current node
                auto [thisID, thisDist] = pq.top();
                pq.pop();

                // if the node was visited before, get out
                // (the shortest distance is already stored in nodeDist)
                if (nodeDist.find(thisID) != nodeDist.end())
                        continue;

                // store the shortest distance to the current node
                nodeDist[thisID] = thisDist;

                // number of nodes exceeded
                if (nodeDist.size() > maxNodes)
                        return;

                SSNode n = getSSNode(thisID);

                // don't go any further if we have reached the maximum length
                int deltaDist = (thisID == srcID) ? 0 : n.getMarginalLength();
                if (thisDist + deltaDist > maxLen)
                        continue;

                // process the right arcs
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++)
                        pq.push(NodeDepth(it->getNodeID(), thisDist + deltaDist));
        }
}

void DBGraph::getSubgraph(NodeID srcID, NodeID dstID, set<NodeID>& nodes)
{
        nodes.clear();

        vector<NodeID> pq;
        pq.push_back(srcID);

        while (!pq.empty()) {
                // get and erase the next node
                NodeID thisID = pq.back();
                pq.pop_back();

                // if the node was visited before, get out
                if (nodes.find(thisID) != nodes.end())
                        continue;

                // store node
                nodes.insert(thisID);

                // target reached
                if (thisID == dstID)
                        continue;

                // process the right arcs
                SSNode n = getSSNode(thisID);
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++)
                        pq.push_back(it->getNodeID());
        }
}

bool DBGraph::findPath(NodeID srcID, NodeID dstID, double minCov, int maxLen) const
{
        set<NodeID> visited;

        priority_queue<NodeDepth, vector<NodeDepth>, NodeDepthComp> pq;
        pq.push(NodeDepth(srcID, 0));

        while (!pq.empty()) {
                // get and erase the current node
                auto [thisID, thisDist] = pq.top();
                pq.pop();

                // if the node was visited before
                if (visited.find(thisID) != visited.end())
                        continue;       // skip node

                // store the shortest distance to the current node
                visited.insert(thisID);

                // target reached?
                if (thisID == dstID)
                        return true;

                SSNode n = getSSNode(thisID);

                // sufficient coverage
                if (n.getAvgCov() <= minCov)
                        continue;

                // don't go any further if we have reached the maximum length
                int deltaDist = (thisID == srcID) ? 0 : n.getMarginalLength();
                if (thisDist + deltaDist > maxLen)
                        continue;

                // process the right arcs
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++)
                        pq.push(NodeDepth(it->getNodeID(), thisDist + deltaDist));
        }

        return false;
}

bool DBGraph::toRome(NodeID srcID, NodeID dstID, int maxLen)
{
        map<NodeID, int> visited;

        priority_queue<NodeDepth, vector<NodeDepth>, NodeDepthComp> pq;
        pq.push(NodeDepth(srcID, 0));

        while (!pq.empty()) {
                // get and erase the current node
                auto [thisID, thisDist] = pq.top();
                pq.pop();

                // if the node was visited before
                if (visited.find(thisID) != visited.end())
                        if (thisID == srcID)
                                return false;   // no backloops allowed
                        else
                                continue;       // skip node

                // store the shortest distance to the current node
                visited[thisID] = thisDist;

                // target reached?
                if (thisID == dstID)
                        continue;

                SSNode n = getSSNode(thisID);

                // don't go any further if we have reached the maximum length
                int deltaDist = (thisID == srcID) ? 0 : n.getMarginalLength();
                if (thisDist + deltaDist > maxLen)
                        return false;

                // if we reach a dead end (different from target), get out
                if (n.numRightArcs() == 0)
                        return false;

                // process the right arcs
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++)
                        pq.push(NodeDepth(it->getNodeID(), thisDist + deltaDist));
        }

        return true;
}

void DBGraph::findTangles(NodeMap<int>& nodeMap,
                          vector<NodePair>& tangles)
{
        tangles.clear();
        int maxLen = 300000, maxNodes = 30;

        // generate candidate tangles, together with the size of the subgraph
        vector<pair<size_t, NodePair>> candT;
        for (NodeID srcID = -getNumNodes(); srcID <= getNumNodes(); srcID++) {
                if ((srcID == 0) || nodeMap[srcID] != 1)
                        continue;
                SSNode n = getSSNode(srcID);
                if (!getSSNode(srcID).isValid() || n.numRightArcs() != 1)
                        continue;

                // get a subgraph from srcID
                map<NodeID, int> nodeDist;
                getSubgraph(srcID, maxLen, maxNodes, nodeDist);

                // filter from this subgraph the unique nodes
                multimap<int, NodeID> uniqueNodes;
                for (auto it : nodeDist) {
                        if (nodeMap[it.first] != 1)
                                continue;
                        if (it.second == 0)     // distance from srcID == 0 ?
                                continue;       // skip direct connections
                        uniqueNodes.insert(make_pair(it.second, it.first));
                }

                // we process nodes from closest to furthest from srcID
                for (auto it : uniqueNodes) {
                        NodeID dstID = it.second;
                        if (!toRome( srcID,  dstID, maxLen) ||
                            !toRome(-dstID, -srcID, maxLen))
                                continue;       // not a tangle
                        NodePair tangle = (srcID < -dstID) ?
                                make_pair(srcID, dstID) :
                                make_pair(-dstID, -srcID);
                        set<NodeID> subgraph;
                        getSubgraph(srcID, dstID, subgraph);
                        candT.push_back(make_pair(subgraph.size(), tangle));
                        break;
                }
        }

        // sort candidate tangles by the length of their subgraph
        sort(candT.begin(), candT.end());
        candT.erase(unique(candT.begin(), candT.end()), candT.end());

        // filter nested tangles (keep the one with the smallest subgraph)
        for (auto el : candT) {
                NodeID srcID = el.second.first;
                NodeID dstID = el.second.second;

                set<NodeID> subgraph;
                getSubgraph(srcID, dstID, subgraph);

                bool overlaps = false;
                for (NodeID nodeID : subgraph) {
                        // only mark internal nodes
                        if ((nodeID == srcID) || (nodeID == dstID))
                                continue;
                        SSNode n = getSSNode(nodeID);
                        if (n.getFlag1()) {
                                overlaps = true;
                                break;
                        }

                        getSSNode(nodeID).setFlag1(true);
                }

                if (!overlaps)
                        tangles.push_back(NodePair(srcID, dstID));
        }

        // reset all flags back to false
        setAllFlags1(false);
}

bool DBGraph::hasUniqueEulerianPath(NodeID srcID, NodeID dstID,
                                    const map<EdgeID, int>& edgeMult) const
{
        // assert that all edge multiplicities are greater than zero
        for (const auto& em : edgeMult)
                assert(em.second > 0);

        // if any node has more than 2 outgoing arcs, there are multiple paths
        // (this is a check that quickly eliminates many candidates)
        map<NodeID, int> numOutArcs;
        for (const auto& em : edgeMult)
                numOutArcs[em.first.first]++;
        for (const auto& na : numOutArcs)
                if (na.second > 2)
                        return false;

        // We use the BEST theorem to count the number of Eulerian *cycles*.
        // Hence, the number of Eulerian paths between srcID and dstID equals
        // the number of Eulerian cyles if we add and edge from srcID to dstID.

        // check if node in-degree == out-degree
        map<NodeID, int> nodeInDegree, nodeOutDegree;
        for (const auto& em : edgeMult) {
                nodeOutDegree[em.first.first] += em.second;
                nodeInDegree[em.first.second] += em.second;
        }
        nodeOutDegree[dstID]++; nodeInDegree[srcID]++;  // add missing link
        if (nodeInDegree != nodeOutDegree)
                return false;                           // necessary condition
        const map<NodeID, int>& nodeDegree = nodeInDegree;

        // construct the (modified) Laplacian matrix Q* obtained by removing an
        // arbitrary row and column (we choose srcID) from the Laplacian matrix
        Matrix<double> Q(nodeDegree.size() - 1, nodeDegree.size() - 1);
        Q.fill(0.0);

        // create a nodeID to index mapping
        map<NodeID, size_t> node2idx; size_t idx = 0;
        for (const auto& nd : nodeDegree)
                if (nd.first != srcID)                   // we removed srcID
                        node2idx[nd.first] = idx++;

        // fill in the elements of Q*
        for (const auto& em : edgeMult) {
                if (em.first.first == em.first.second)
                        continue;       // skip self-loops
                if (em.first.second == srcID)
                        continue;
                size_t j = node2idx[em.first.second];
                Q(j, j) += em.second;   // indegree without self-loops
                if (em.first.first == srcID)
                        continue;
                size_t i = node2idx[em.first.first];
                Q(i, j) = -em.second;   // minus number of edges from i to j

        }

        /*cout << " ================ " << endl;
        cout << "FROM: " << srcID << " -- TO: " << dstID << endl;
        for (size_t i = 0; i < nodeDegree.size() - 1; i++) {
                for (size_t j = 0; j < nodeDegree.size() - 1; j++)
                        cout << Q(i, j) << " ";
                cout << "\n";
        }
        cout << "Determinant: " << Q.determinant() << endl;*/

        // compute the number of Eulerian paths using the BEST theorem:
        // nEurCycle = det(Q) * PROD_i (n_i - 1)! / PROD_j (e_j!)

        // compute numerator and denominator separately
        size_t num = round(Q.determinant()); size_t den = 1;

        // get the values for the numerator
        multiset<int> numValSet;
        for (const auto& nd : nodeDegree)
                numValSet.insert(nd.second - 1);

        // get the values for the denominator
        multiset<int> denValSet;
        for (auto el : edgeMult)
                denValSet.insert(el.second);

        // compute numerator and denominator by grouping factorial values such
        // that overflows are avoided
        while (!(numValSet.empty() && denValSet.empty())) {
                int a = numValSet.empty() ? 1 : *numValSet.rbegin();
                int b = denValSet.empty() ? 1 : *denValSet.rbegin();
                if (a > b)
                        num *= Util::factRatio(a, b);
                else
                        den *= Util::factRatio(b, a);
                if (!numValSet.empty()) // remove largest element
                        numValSet.erase(prev(numValSet.end()));
                if (!denValSet.empty()) // remove largest element
                        denValSet.erase(prev(denValSet.end()));
        }

        //cout << "Number of Eulerian paths: " << (double)num/double(den) << endl;
        return (num == den);
}

vector<NodeID> DBGraph::getEulerianPath(NodeID srcID,
                                        const map<EdgeID, int>& edgeMult) const
{
        // run Hierholzer's algorithm
        vector<NodeID> stack, path;
        stack.push_back(srcID);
        map<EdgeID, int> edgeMultWork = edgeMult;

        while (!stack.empty()) {
                NodeID curr = stack.back();

                // try and extend the current path further
                NodeID next = 0;
                SSNode n = getSSNode(curr);
                for (ArcIt a = n.rightBegin(); a != n.rightEnd(); a++) {
                        EdgeID e(curr, a->getNodeID());
                        auto it = edgeMultWork.find(e);
                        if ( (it == edgeMultWork.end()) || (it->second == 0) )
                                continue;       // edge exhausted
                        next = a->getNodeID();
                        break;
                }

                if (next != 0) {        // extension successful
                        edgeMultWork[EdgeID(curr, next)]--;
                        stack.push_back(next);
                } else {                // backtracking
                        path.push_back(curr);
                        stack.pop_back();
                }
        }

        reverse(path.begin(), path.end());
        return path;
}

NodeChain DBGraph::resolveTangle(const NodePair& tangle,
                                 BiNodeMap<PathInfo>& pathInfo,
                                 NodeMap<int>& nodeMult, EdgeMap<int>& edgeMult,
                                 double MCD)
{
        NodeID srcID = tangle.first;
        NodeID dstID = tangle.second;

        // get all nodes
        set<NodeID> sgNodes;
        getSubgraph(srcID, dstID, sgNodes);

        // populate sgEdgeMult
        map<EdgeID, int> sgEdgeMult;
        for (NodeID id : sgNodes) {
                SSNode n = getSSNode(id);
                if (id == srcID)
                        continue;       // skip incoming arcs of srcID
                if (nodeMult[id] == 0)
                        continue;       // skip zero-mult nodes/edges

                // handle all incoming arcs (this suffices)
                for (ArcIt a = n.leftBegin(); a != n.leftEnd(); a++) {
                        EdgeID e(a->getNodeID(), id);
                        if (edgeMult[e] > 0)
                                sgEdgeMult[e] = edgeMult[e];
                }
        }

        //cout << "Solving tangle from: " << srcID << " to " << dstID << endl;

        // is the Euler path unique?
        if (hasUniqueEulerianPath(srcID, dstID, sgEdgeMult)) {
                NodeChain r(getEulerianPath(srcID, sgEdgeMult));
                //cout << "Reduction via UEP: " << r << endl;
                return r;
        }

        // take into account the source and destination chain
        NodeChain srcChain = pathInfo[srcID].getMostLikelyPath(MCD).toNodeChain();
        if (srcChain.empty())   // prevent pathological case
                srcChain.push_back(srcID);
        NodeChain dstChain = pathInfo[-dstID].getMostLikelyPath(MCD).toNodeChain();
        dstChain.revCompl();
        if (dstChain.empty())   // prevent pathological case
                dstChain.push_back(dstID);

        /*for (auto el : sgEdgeMult)
                cout << el.first.first << " --> " << el.first.second
                     << "\t(" << el.second << ")\n";
        cout << "srcChain: " << srcChain << endl;
        cout << "dstChain: " << dstChain << endl;*/

        // process source chain
        for (size_t i = 1; i < srcChain.size(); i++) {
                EdgeID eID(srcChain[i-1], srcChain[i]);
                auto it = sgEdgeMult.find(eID);
                if (it != sgEdgeMult.end()) {
                        it->second--;
                        if (it->second == 0)
                                sgEdgeMult.erase(it);
                } else {
                        // multiplicities are not correctly infered
                        return NodeChain();
                }

                if (sgEdgeMult.empty()) { // tangle fully consumed
                        NodeChain r(srcChain.begin(), srcChain.begin() + i + 1);
                        //cout << "Reduction via srcChain: " << r << endl;
                        return r;
                }
        }

        // process destination chain
        for (size_t i = dstChain.size(); i--> 1; ) {
                EdgeID eID(dstChain[i-1], dstChain[i]);
                auto it = sgEdgeMult.find(eID);
                if (it != sgEdgeMult.end()) {
                        it->second--;
                        if (it->second == 0)
                                sgEdgeMult.erase(it);
                } else {
                        // multiplicities are not correctly infered
                        return NodeChain();
                }

                if (sgEdgeMult.empty()) { // tangle fully consumed
                        NodeChain r = srcChain;
                        while (srcChain[srcChain.size()-1] != dstChain[i-1]) { // bugfix when flow might not hold despite multiplicity estimation
                                i++; // move back up the dest chain until we can attach to src chain
                        }
                        r.insert(r.end(), dstChain.begin() + i, dstChain.end());
                        //cout << "Reduction via src + dst: " << r << endl;
                        return r;
                }
        }

        /*for (auto el : sgEdgeMult)
                cout << el.first.first << " --> " << el.first.second
                     << "\t(" << el.second << ")\n";*/

        if (hasUniqueEulerianPath(srcChain.back(), dstChain.front(), sgEdgeMult)) {
                NodeChain p(getEulerianPath(srcChain.back(), sgEdgeMult));
                NodeChain r = srcChain;
                r.pop_back();
                r.insert(r.end(), p.begin(), p.end());
                r.pop_back();
                r.insert(r.end(), dstChain.begin(), dstChain.end());

                //cout << "Reduction via partial UEP: " << r << endl;
                return r;
        }

        return NodeChain();
}

void DBGraph::resolveTangles(LibraryContainer& library,
                             NodeMap<int>& nodeMult,
                             EdgeMap<int>& edgeMult,
                             vector<NodePair>& tangles,
                             std::vector<NodeChain>& reductions,
                             NodeMap<int>& redNodeMult,
                             EdgeMap<int>& redEdgeMult)
{
        const float MCD = 0.25;         // FIXME! Parameter!
        const int minCount = 2;         // FIXME! Parameter!

        reductions.clear();
        redNodeMult.clear();
        redEdgeMult.clear();

        BiNodeMap<bool> noi;            // indicate nodes of interest
        for (const auto& it : tangles) {
                noi[ it.first]  = true;
                noi[-it.second] = true;
        }

        // build the path info from the aligned reads
        BiNodeMap<PathInfo> pathInfo;
        addReadInfo(library, noi, pathInfo);
        filterPathInfo(pathInfo, MCD, minCount);

        NodeChainSet red;
        for (auto tangle : tangles) {
                NodeChain r = resolveTangle(tangle, pathInfo,
                                            nodeMult, edgeMult, MCD);
                if (!r.empty()) {
                        red.insert(r);
                        r.revCompl();
                        red.insert(r);
                }
        }

        // aggregate the reductions into larger contigs
        BiNodeMap<Contig> contigs;
        BiNodeMap<int> node2contig;

        buildContigs(red, contigs, node2contig);
        cout << "Clustered reductions into " << contigs.size()/2 << " contigs\n";

        reductions.clear();
        for (auto& [ignore, cnc] : contigs) {
                NodeChain nc = cnc.toNodeChain();
                if (nc.size() <= 1)
                        continue;
                if (nc < nc.getRevCompl())
                        reductions.emplace_back(nc);
        }

        // set the node and edge multiplicities for the internal nodes
        for (const auto& r : reductions) {
                if (r.size() < 3)
                        continue;
                for (size_t i = 1; i < r.size() - 1; i++)
                        redNodeMult[NodeRep(r[i])]++;
                for (size_t i = 1; i < r.size() - 2; i++)
                        redEdgeMult[EdgeRep(r[i], r[i+1])]++;
        }
}
