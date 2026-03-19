/******************************************************************************
 *   Copyright (C) 2022 Jan Fostier (jan.fostier@ugent.be)                    *
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

#include "nodetracker.h"
#include "dbgraph.h"

using namespace std;

NodeLength NodeTracker::getNodeLength(NodeID oldID) const
{
        auto it = log.find(oldID);
        return (it == log.end()) ?
                dBG.getSSNode(oldID).getMarginalLength() :
                it->second.margLen;
}

void NodeTracker::build(const vector<NodeChain>& reductions,
                        const NodeMap<int>& nodeMult,
                        const NodeMap<int>& redNodeMult)
{
        log.clear();

        // build the nodeTracker infolog
        for (const auto& r : reductions) {
                // compute the node offsets
                NodeLength offset = 0;
                for (size_t i = 0; i < r.size(); i++) {
                        const int enm = get(nodeMult, r[i]);
                        const int rnm = get(redNodeMult, r[i]);
                        NodeLength mnl = dBG.getSSNode(r[i]).getMarginalLength();
                        if (enm == 1) {
                                assert(rnm == 1);
                                // unique node becomes part of 'bigger' node
                                logEvent(r[i], r.front(), offset, mnl);
                        } else if (rnm == enm) {
                                // repeated node will be deleted
                                logEvent(r[i], 0, 0, mnl);
                        }
                        offset += mnl;
                }

                offset = 0;
                for (size_t i = r.size(); i--> 0; ) {
                        const int enm = get(nodeMult, -r[i]);
                        const int rnm = get(redNodeMult, -r[i]);
                        NodeLength mnl = dBG.getSSNode(-r[i]).getMarginalLength();
                        if (enm == 1) {
                                assert(rnm == 1);
                                // unique node becomes part of 'bigger' node
                                logEvent(-r[i], -r.front(), offset, mnl);
                        } else if (rnm == enm) {
                                // repeated node will be deleted
                                logEvent(-r[i], 0, 0, mnl);
                        }
                        offset += mnl;
                }
        }
}

void NodeTracker::build(const vector<NodeChain>& reductions)
{
        log.clear();

        // build the nodeTracker infolog
        for (const auto& r : reductions) {
                // compute the node offsets
                NodeLength offset = 0;
                for (size_t i = 0; i < r.size(); i++) {
                        NodeLength mnl = dBG.getSSNode(r[i]).getMarginalLength();
                        // node becomes part of 'bigger' node
                        logEvent(r[i], r.front(), offset, mnl);
                        offset += mnl;
                }

                offset = 0;
                for (size_t i = r.size(); i--> 0; ) {
                        NodeLength mnl = dBG.getSSNode(-r[i]).getMarginalLength();
                        // node becomes part of 'bigger' node
                        logEvent(-r[i], -r.front(), offset, mnl);
                        offset += mnl;
                }
        }
}

void NodeTracker::updateNodeChain(NodeChain& oldC, int& begAln, int& endAln) const
{
        if (log.empty())
                return;

        // create a new node chain with updated nodeIDs
        NodeChain newC(oldC);
        vector<NodeLength> offs(oldC.size());
        vector<int> mlen(oldC.size());
        vector<bool> anch(oldC.size());
        size_t la = oldC.size(), ra = 0;  // leftmost and rightmost anchor node

        for (size_t i = 0; i < newC.size(); i++) {
                newC[i] = getPresentNodeID(oldC[i]);
                offs[i] = getNodeOffset(oldC[i]);
                mlen[i] = getNodeLength(oldC[i]);
                anch[i] = isAnchorNode(oldC[i]);
                if (anch[i]) {
                        la = min(la, i);
                        ra = max(ra, i);
                }
        }

        //cout << "ORIG: " << oldC << "  B: " << begAln << "  E: " << endAln << endl;
        //cout << "TRAN: " << newC << endl;

        // a) check if the tail needs to be removed
        if (la != oldC.size()) {
                // we know that newC[ra] is the rightmost anchor node
                size_t j = ra;
                for ( ; j < newC.size() - 1; j++)
                        if (!dBG.edgeExists(EdgeRep(newC[j], newC[j+1])))
                                break;       // invalid path found

                // nodes [j+1, end[ will be deleted
                int lenToDelete = 0;
                for (size_t k = j; k < newC.size() - 1; k++)
                        lenToDelete += mlen[k];
                lenToDelete += endAln;

                endAln = offs[j] + lenToDelete;

                if ((j + 1) < newC.size())
                        newC.erase(newC.begin() + j + 1, newC.end());
        }

        //cout << "TAIL: " << newC << "  E: " << endAln << endl;

        // b) check if the head needs to be removed
        if (la != oldC.size()) {
                // at this point, newC[la] is the leftmost anchor node
                size_t j = la;
                for ( ; j > 0; j--)
                        if (!dBG.edgeExists(EdgeRep(newC[j-1], newC[j])))
                                break;       // invalid path found

                // nodes [0, j[ will be deleted
                int lenToDelete = 0;
                for (size_t k = 0; k < j; k++)
                        lenToDelete += mlen[k];
                lenToDelete -= begAln;

                begAln = offs[j] - lenToDelete;

                if (j > 0) {    // delete nodes [0, j[
                        newC.erase(newC.begin(), newC.begin() + j);
                        anch.erase(anch.begin(), anch.begin() + j);
                }
        }

        //cout << "HEAD: " << newC << "  B: " << begAln << endl;
        //if (newC.front() == newC.back())
        //        cout << "E-B: " << endAln - begAln << endl;

        // c) remove all nodes in between successive, identical anchor nodes
        for (size_t i = newC.size(), j = 0; j < newC.size(); j++) {
                if (!anch[j])
                        continue;       // not an anchor node

                // remove everything in between identical anchor nodes
                if ((i != newC.size()) && (newC[i] == newC[j])) {
                        newC.erase(newC.begin() + i, newC.begin() + j);
                        anch.erase(anch.begin() + i, anch.begin() + j);
                        j = i;
                        continue;
                }

                i = j;  // open interval [i, ...]
        }

        //cout << "FINAL:" << newC << endl;
        //int tmp; cin >> tmp;

        if (!dBG.pathExists(newC))
                newC.clear();
        oldC = newC;
}
