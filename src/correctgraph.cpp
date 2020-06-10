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

#include "dbgraph.h"

using namespace std;

void DBGraph::removeNode(NodeID nodeID)
{
        SSNode node = getSSNode(nodeID);
        if (!node.isValid())            // node already removed
                return;

        for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++) {
                SSNode leftNode = getSSNode(it->getNodeID());
                if (abs(leftNode.getNodeID()) == abs(node.getNodeID()))
                        continue;       // skip self-loops or palindromic arcs
                leftNode.deleteRightArc(node.getNodeID());
                numValidArcs--;
        }

        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                SSNode rightNode = getSSNode ( it->getNodeID() );
                if (abs(rightNode.getNodeID()) == abs(node.getNodeID()))
                        continue;       // skip self-loops or palindromic arcs
                rightNode.deleteLeftArc(node.getNodeID());
                numValidArcs--;
        }

        numValidArcs -= node.numRightArcs();
        node.deleteAllRightArcs();
        numValidArcs -= node.numLeftArcs();
        node.deleteAllLeftArcs();

        numValidNodes--;
        node.invalidate();
}

void DBGraph::removeArc(NodeID leftID, NodeID rightID)
{
        if (getSSNode(leftID).deleteRightArc(rightID))
                numValidArcs--;

        if (getSSNode(rightID).deleteLeftArc(leftID))
                numValidArcs--;
}

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

        cout << "\tRemoved " << initNumValidNodes - numValidNodes << " nodes\n";
        cout << "\tRemoved " << initNumValidArcs - numValidArcs << " arcs\n";
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

void DBGraph::convertNodesToString(const vector<NodeID> &nodeSeq,
                                   string &output)
{
        output.clear();

        if (nodeSeq.empty())
                return;

        size_t size = 0;
        for (size_t i = 0; i < nodeSeq.size(); i++)
                size += getSSNode(nodeSeq[i]).length();
        size -= (nodeSeq.size() - 1) * (Kmer::getK() - 1);

        output = getSSNode(nodeSeq[0]).getSequence();
        for (size_t i = 1; i < nodeSeq.size(); i++)
                output.append(getSSNode(nodeSeq[i]).getSequence().substr(Kmer::getK() - 1));

        assert(size == output.size());
}

void DBGraph::concatenateAroundNode(NodeID seedID, vector<NodeID>& nodeListv)
{
        nodeListv.clear();

        SSNode seed = getSSNode(seedID);
        if (!seed.isValid())
                return;

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

        string str;
        convertNodesToString(nodeListv, str);

        Coverage newCov = front.getCov();
        for (size_t i = 1; i < nodeListq.size(); i++) {
                SSNode n = getSSNode(nodeListq[i]);

                newCov += n.getCov();
                n.deleteAllLeftArcs();
                n.deleteAllRightArcs();
                n.invalidate();
        }

        front.setCov(newCov);
        front.setSequence(str);
        numValidNodes -= nodeListq.size() - 1;
        numValidArcs -= 2*(nodeListq.size() - 1);
}

bool DBGraph::concatenateNodes()
{
        size_t numConcatenations = 0;

        for (NodeID seedID = 1; seedID <= numNodes; seedID++) {
                vector<NodeID> concatenation;
                concatenateAroundNode(seedID, concatenation);

                if (!concatenation.empty())
                        numConcatenations += concatenation.size() - 1;
        }

        cout << "\tConcatenated " << numConcatenations << " nodes" << endl;

        size_t countTotal = 0;
        for (NodeID seedID = 1; seedID <= numNodes; seedID++)
                if (getSSNode(seedID).isValid())
                        countTotal++;

        return (numConcatenations > 0);
}

void DBGraph::removeNodes(const std::vector<NodeRep>& nodes)
{
        for (const auto& el : nodes)
                removeNode(el);
}

void DBGraph::removeEdges(const std::vector<EdgeRep>& edges)
{
        for (const auto& el : edges)
                removeArc(el.getSrcID(), el.getDstID());
}
