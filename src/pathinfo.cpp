/******************************************************************************
 *   Copyright (C) 2014 - 2023 Jan Fostier (jan.fostier@ugent.be)             *
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

#include <algorithm>
#include "nodechain.h"
#include "pathinfo.h"

using namespace std;

// ============================================================================
// AUXILIARY ROUTINES
// ============================================================================

bool sortSecondReverse(const NodeCount &a, const NodeCount &b) {
        return (a.second > b.second);
}

// ============================================================================
// CHAIN COLLECTION
// ============================================================================

void ChainColl::addChain(const NodeChain& chain)
{
        // add the node chain
        if (chain.empty())
                return;

        set<NodeChain>::iterator it = find(chain);
        if (it == end())
                insert(chain);
        else
                const_cast<NodeChain&>(*it).incrCount(chain.getCount());
}

// ============================================================================
// PATHINFO (private)
// ============================================================================

PathInfo::PathInfo(const ChainColl& cc)
{
        // build direct consensus node chains
        for (const NodeChain& nc : cc)
                if (dirCNC.empty() || !dirCNC.back().mergeNodeChain(nc))
                        dirCNC.push_back(nc);
}

BiNodeMap<int> PathInfo::getNumChainsPerNode() const
{
        BiNodeMap<int> result;
        for (const ConsNodeChain& dir : dirCNC) {
                // make the counts per consensus chain
                BiNodeMap<int> cpn;
                for (const auto& [id, count] : dir)
                        cpn[id] = max<int>(cpn[id], count);
                // aggregate into global count
                for (const auto& [id, count] : cpn)
                        result[id] += count;
        }

        return result;
}

set<NodeID> PathInfo::filterMCD(const vector<NodeCount>& nodeCounts, float MCD)
{
        // sort the counts high-to-low
        vector<NodeCount> countsSorted = nodeCounts;
        std::sort(countsSorted.begin(), countsSorted.end(), sortSecondReverse);

        // nodes are filtered from index [i+1, end[ when the cumulative count
        // of nodes [i+1, end[ has a CD <= MCD w.r.t. the count of node[i]

        // compute the cumulative count cumCount[i] = sum(count[i+1, end[)
        vector<int> cumCount(countsSorted.size(), 0);
        int sum = 0;
        for (size_t i = countsSorted.size(); i--> 0; ) {
                cumCount[i] = sum;
                sum += countsSorted[i].second;
        }

        // filter the nodes
        set<NodeID> result;
        for (size_t i = 0; i < countsSorted.size(); i++) {
                if (CD(countsSorted[i].second, cumCount[i]) <= MCD) {
                        for (size_t j = i + 1; j < countsSorted.size(); j++)
                                result.insert(countsSorted[j].first);
                        return result;
                }
        }

        return result;
}

// ============================================================================
// PATHINFO (public)
// ============================================================================

void PathInfo::filterDirPaths(float MCD)
{
        vector<tuple<size_t, size_t, size_t>> stack(1, {0, 0, dirCNC.size()});
        while (!stack.empty()) {
                auto [d, b, e] = stack.back();
                stack.pop_back();

                //cout << "d = " << d << ", [" << b << ", " << e << "]" << endl;

                vector<NodeCount> nodeCounts;
                for (size_t i = b; i < e; i++) {
                        if (d >= dirCNC[i].size())
                                continue;
                        auto [id, count] = dirCNC[i][d];
                        if (nodeCounts.empty() || nodeCounts.back().first != id)
                                nodeCounts.emplace_back(id, count);
                        else
                                nodeCounts.back().second += count;
                }

                set<NodeID> toFilter = filterMCD(nodeCounts, MCD);

                vector<tuple<size_t, size_t, size_t>> next;
                NodeID prevID = 0;
                for (size_t i = b; i < e; i++) {
                        if (d >= dirCNC[i].size())
                                continue;
                        auto [id, count] = dirCNC[i][d];
                        if (toFilter.find(id) != toFilter.end()) {
                                dirCNC[i] = ConsNodeChain(dirCNC[i].begin(),
                                                          dirCNC[i].begin()+d);
                                continue;
                        }

                        if (next.empty() || prevID != id) {
                                next.emplace_back(d+1, i, i+1);
                                prevID = id;
                        } else
                                get<2>(next.back())++;

                }

                reverse(next.begin(), next.end());
                for (const auto e : next)
                        stack.push_back(e);
        }

        // compress dirCNC
        vector<ConsNodeChain> copy;
        swap(copy, dirCNC);
        for (const auto& nc : copy) {
                if (dirCNC.empty() || !dirCNC.back().mergeNodeChain(nc))
                        dirCNC.push_back(nc);
        }
}

BiNodeMap<int> PathInfo::getUniqueDirDst(const BiNodeMap<bool>& unique) const
{
        BiNodeMap<int> result;
        for (const ConsNodeChain& dir : dirCNC) {
                // make the counts per consensus chain
                BiNodeMap<int> cpn;
                for (size_t i = 1; i < dir.size(); i++) {
                        const auto& [id, count] = dir[i];
                        if (get(unique, -id))
                                cpn[id] = max<int>(cpn[id], count);
                }
                // aggregate into global count
                for (const auto& [id, count] : cpn)
                        result[id] += count;
        }

        return result;
}

BiNodeMap<int> PathInfo::getUniquePerDst(const BiNodeMap<bool>& unique) const
{
        BiNodeMap<int> result;
        for (const auto [readID, p] : per) {
                // find the destinations per paired-end read
                set<NodeID> dst;
                for (const auto id : p.first)
                        if (get(unique, -id))
                                dst.insert(id);
                for (const auto id : p.second)
                        if (get(unique, -id))
                                dst.insert(id);
                // aggregate into global count
                for (const auto& id : dst)
                        result[id]++;
        }

        return result;
}

BiNodeMap<ReadIDSet> PathInfo::getDstReadSet(const BiNodeMap<bool>& noi) const
{
        // collect, per dstID, a set of readIDs that map to it
        BiNodeMap<ReadIDSet> result;
        for (const auto& [readID, readPair] : per) {
                for (NodeID id : readPair.first)
                        if (get(noi, id))
                                result[id].insert(readID);
                for (NodeID id : readPair.second)
                        if (get(noi, id))
                                result[id].insert(readID);
        }

        return result;
}

ConsNodeChain PathInfo::getMostLikelyPath(float MCD) const
{
        ConsNodeChain result;

        for (size_t i = 0; ; i++) {
                // make a map of node counts at offset i
                BiNodeMap<int> nodeCounts;
                for (const ConsNodeChain& nc : dirCNC) {
                        if (i < nc.size()) {
                                auto [id, count] = nc[i];
                                nodeCounts[id] += count;
                        }
                }

                // all paths have reached the end
                if (nodeCounts.empty())
                        return result;

                // compute conflict degree
                int sumVal = 0, maxVal = -1; NodeID maxKey = 0;
                for (const auto& [id, count] : nodeCounts) {
                        sumVal += count;
                        if (count > maxVal) {
                                maxVal = count;
                                maxKey = id;
                        }
                }

                // check maximum conflict degree
                if (CD(maxVal, sumVal - maxVal) > MCD)
                        return result;

                result.emplace_back(maxKey, maxVal);
        }
}

bool PathInfo::isCompatible(const NodeChain& ncoi, float MCD) const
{
        for (size_t i = 0; i < ncoi.size(); i++) {
                // count the number of corresponding and discordant chains
                int nCorr = 0, nDisc = 0;
                for (const ConsNodeChain& nc : dirCNC) {
                        if (i < nc.size()) {
                                if (ncoi[i] == nc[i].first)
                                        nCorr += nc[i].second;
                                else
                                        nDisc += nc[i].second;
                        }
                }

                // all paths have reached the end
                if (nCorr == 0 && nDisc == 0)
                        return true;

                // check maximum conflict degree
                if (CD(nCorr, nDisc) > MCD)
                        return false;
        }

        return true;
}

bool PathInfo::isUnique(float MCD) const
{
        // This function checks the uniqueness of a node based on its reads.
        // We first count how many reads N pass through this node. If more than
        // MCD * N reads diverge from the main path, uniqueness is rejected.
        int numChains = 0;
        for (const ConsNodeChain& nc : dirCNC)
                if (!nc.empty())
                        numChains += nc.front().second;

        for (size_t i = 1; ; i++) {
                // make a map of node counts at offset i
                BiNodeMap<int> nodeCounts;
                for (const ConsNodeChain& nc : dirCNC) {
                        if (i < nc.size()) {
                                auto [id, count] = nc[i];
                                nodeCounts[id] += count;
                        }
                }

                // all paths have reached the end
                if (nodeCounts.empty())
                        return true;

                // compute conflict degree
                int sumVal = 0, maxVal = -1;
                for (const auto [id, count] : nodeCounts) {
                        sumVal += count;
                        maxVal = max<int>(maxVal, count);
                }

                // check maximum conflict degree (relative to numChains!)
                if (CD(numChains, sumVal - maxVal) > MCD)
                        return false;
        }

        return true;
}

void PathInfo::addUniqueReads(std::vector<NodeChain>& reads,
                              std::set<size_t>& readIDs, bool revCompl) const
{
        for (const auto& [readID, pe] : per) {
                if (!readIDs.insert(readID).second)
                        continue;       // paired-end read already added

                const NodeChain& L = pe.first;
                if (!L.empty())
                        reads.emplace_back(revCompl ? L.getRevCompl() : L);

                const NodeChain& R = pe.second;
                if (!R.empty())
                        reads.emplace_back(revCompl ? R.getRevCompl() : R);
        }
}

void PathInfo::buildPathFinder(BiNodeMap<PathFinder>& pf,
                               set<size_t>& readIDs) const
{
        for (const auto& [readID, pe] : per) {
                if (!readIDs.insert(readID).second)
                        continue;       // paired-end read already in graph

                const NodeChain& L = pe.first;
                for (size_t i = 1; i + 1 < L.size(); i++) {
                        cout << i << endl;
                        if (pf.find(L[i]) == pf.end())
                                continue;       // not of interest
                        pf[L[i]].addChain(NodeChain(L.begin(), L.begin() + i), L[i+1]);
                        cout << "Adding" << endl;
                }

                const NodeChain& R = pe.second;
                for (size_t i = 1; i + 1 < R.size(); i++) {
                        if (pf.find(R[i]) == pf.end())
                                continue;       // not of interest
                        pf[R[i]].addChain(NodeChain(R.begin(), R.begin() + i), R[i+1]);
                }
        }
}
