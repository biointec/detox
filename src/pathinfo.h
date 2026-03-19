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

#ifndef PATHINFO_H
#define PATHINFO_H

#include <vector>
#include <set>

#include "ssnode.h"

// ============================================================================
// DEFINITIONS
// ============================================================================

typedef std::pair<NodeID, int> NodeCount;
typedef std::set<size_t> ReadIDSet;

class PathFinder;

// ============================================================================
// CHAIN COLLECTION
// ============================================================================

class ChainColl : public std::set<NodeChain>
{
public:
        /**
         * Add a node chain to the collection
         * @param chain Chain to add
         */
        void addChain(const NodeChain& chain);
};

// ============================================================================
// PATH INFORMATION CLASS
// ============================================================================

class PathInfo
{
//private:
        public:
        std::vector<ConsNodeChain> dirCNC;      // direct consensus node chain
        std::vector<std::pair<size_t, std::pair<NodeChain, NodeChain>>> per;

        float CD(int nCorr, int nDisc) const {
                return (nCorr == 0) ? 1.0f : float(nDisc) / float(nCorr+nDisc);
        }

        /**
         * Get the number of **different** chains that run through a node
         * @return A map containing <NodeID, int> entries
         */
        BiNodeMap<int> getNumChainsPerNode() const;

        /**
         * Given a vector of node counts and a maximum conflict degree,
         * return a vector of nodes to filter (= minority nodes)
         * @param nodeCounts Vector of node counts
         * @param MCD Maximum conflict degree
         * @return A set of nodes to filter
         */
        std::set<NodeID> filterMCD(const std::vector<NodeCount>& nodeCounts,
                                   float MCD);

public:
        PathInfo() {}

        PathInfo(const ChainColl& cc);

        void addPer(const NodeChain& L, const NodeChain& R, size_t readID) {
                per.emplace_back(make_pair(readID, make_pair(L, R)));
        }

        /**
         * Filter minority paths from the consensus direct node chains
         * @param MCD Maximum conflict degree
         */
        void filterDirPaths(float MCD);

        /**
         * Get unique direct target nodes and the number of chains to each node
         * @param unique Node map that flags unique nodes
         * @return BiNodeMap with unique direct targets
         */
        BiNodeMap<int> getUniqueDirDst(const BiNodeMap<bool>& unique) const;

        /**
         * Get unique PER target nodes and the number of chains to each node
         * @param unique Node map that flags unique nodes
         * @return BiNodeMap with unique PER targets
         */
        BiNodeMap<int> getUniquePerDst(const BiNodeMap<bool>& unique) const;

        /**
         * Get a set of readIDs for each DIR/PER destination node
         * @param noi Node map that flags destination nodes
         * @return BiNodeMap with a set of readIDs per DIR/PER destination node
         */
        BiNodeMap<std::set<size_t>>
        getDstReadSet(const BiNodeMap<bool>& noi) const;

        /**
         * Add the paired-end reads to a reads container
         * @param reads Reads (input/output)
         * @param readIDs Set of read identifiers in graph (input/output)
         * @param revCompl True if the reads should be reverse-complemented
         */
        void addUniqueReads(std::vector<NodeChain>& reads,
                            std::set<size_t>& readIDs,
                            bool revCompl = false) const;

        /**
         * @brief Build path finder information
         * @param pf Path finder data structure
         * @param readIDs Set of read identifiers in graph (input/output)
         */
        void buildPathFinder(BiNodeMap<PathFinder>& pf,
                             std::set<size_t>& readIDs) const;

        /**
         * Get the longest, most probable path that does not exceed MCD
         * @param MCD The maximum conflict degree
         * @return A consensus node chain
         */
        ConsNodeChain getMostLikelyPath(float MCD) const;

        /**
         * Check whether a given path does not exceed the MCD
         * @param nc Node Chain
         * @param MCD The maximum conflict degree
         * @return True of false
         */
        bool isCompatible(const NodeChain& nc, float MCD) const;

        /**
         * Check whether a node is unique based on path information
         * @param MCD Maximum conflict degree
         * @return True of false
         */
        bool isUnique(float MCD) const;

        /**
         * @brief Print the direct and paired-end paths
         */
        void printInfo() const {
                for (auto it : dirCNC)
                        std::cout << it << std::endl;
        }

        /**
         * @brief Print the direct and paired-end paths (reverse complemented)
         */
        void printReverseInfo() const {
                for (auto it : dirCNC)
                        std::cout << it.getRevCompl() << std::endl;
        }

        /**
         * Remove all data
         */
        void clear() {
                dirCNC.clear();
                per.clear();
        }
};

// ============================================================================
// TARGET INFORMATION CLASS
// ============================================================================

class NodeChainSet : public std::set<NodeChain>
{
private:

public:
        // keep the base class constructor
        using std::set<NodeChain>::set;

        // keep base class implementations of begin() and end()
        using std::set<NodeChain>::begin;
        using std::set<NodeChain>::end;
        using std::set<NodeChain>::insert;

        /**
         * Get an iterator to the first node chain that starts with id
         * @param id Node identifier
         * @return An iterator to the first node chain that starts with id
         */
        std::set<NodeChain>::iterator begin(NodeID id) {
                return lower_bound(NodeChain({id}));
        }

        /**
         * Get an iterator to the first node chain that starts with id
         * @param id Node identifier
         * @return An iterator to the first node chain that starts with id
         */
        std::set<NodeChain>::const_iterator begin(NodeID id) const {
                return lower_bound(NodeChain({id}));
        }

        /**
         * Get an iterator past the last node chain that starts with id
         * @param id Node identifier
         * @return An iterator past the last node chain that starts with id
         */
        std::set<NodeChain>::iterator end(NodeID id) {
                return lower_bound(NodeChain({id + 1 }));
        }

        /**
         * Get an iterator past the last node chain that starts with id
         * @param id Node identifier
         * @return An iterator past the last node chain that starts with id
         */
        std::set<NodeChain>::const_iterator end(NodeID id) const {
                return lower_bound(NodeChain({id + 1 }));
        }

        /**
         * @brief Insert a nodeChain, update counter if element already exists
         * @param el Element to insert
         * @return std::pair<iterator,bool>
         */
        std::pair<iterator,bool> insert (const value_type& el)
        {
                auto [it, flag] = std::set<NodeChain>::insert(el);
                int maxCount = std::max<int>(it->getCount(), el.getCount());
                const_cast<NodeChain&>(*it).setCount(maxCount);
                return make_pair(it, flag);
        }
};

// ============================================================================
// CONTIG CLUSTER CLASS
// ============================================================================

class Contig : public ConsNodeChain {

private:


public:
        Contig() {};
        Contig(const ConsNodeChain& nc) : ConsNodeChain(nc) {}

        /**
         * Operator<< overloading
         * @param out Output file stream (input)
         * @param ncc ContigCluster to display
         * @return Output file stream
         */
        /*friend std::ostream &operator<<(std::ostream &out,
                                        const ContigCluster& nc);*/
};

// ============================================================================
// PATHFINDER CLASS
// ============================================================================

class PathFinder {

private:
        ChainColl rawData;

        std::vector<ConsNodeChain> upstrm;      // upstream paths
        std::vector<NodeID> next;               // next node

        float CD(int nCorr, int nDisc) const {
                return (nCorr == 0) ? 1.0f : float(nDisc) / float(nCorr+nDisc);
        }

public:
        PathFinder() {
                /*upstrm.emplace_back(ConsNodeChain({{2,10}, {1,1}}));
                upstrm.emplace_back(ConsNodeChain({{3,1}, {1,4}}));
                upstrm.emplace_back(ConsNodeChain({{4,1}}));

                for (size_t i = 0; i < upstrm.size(); i++)
                        std::reverse(upstrm[i].begin(), upstrm[i].end());

                next.emplace_back(-1);
                next.emplace_back(-3);
                next.emplace_back(-2);*/
        }

        void addChain(const NodeChain& upstream, NodeID nextID)
        {
                NodeChain copy = upstream;
                copy.reverse();
                copy.push_back(nextID);
                rawData.addChain(copy);
        }

        void buildIndex() {
                // build direct consensus node chains
                for (const NodeChain& nc : rawData) {
                        NodeChain copy = nc;
                        NodeID nextID = copy.back();
                        copy.pop_back();

                        if (upstrm.empty() || next.back() != nextID || !upstrm.back().mergeNodeChain(copy)) {
                                upstrm.emplace_back(copy);
                                next.push_back(nextID);
                        }
                }
        }

        /**
         * Predict the next node given a path of previously visited nodes
         * @param nc A chain of previously visited nodes
         * @param MCD Maximum conflict degree
         * @return Node identifier to visit next
         */
        NodeID getNextID(const NodeChain& nc, float MCD)
        {
                NodeChain path = nc;
                path.pop_back();
                std::reverse(path.begin(), path.end());

                std::vector<ConsNodeChain::iterator> it(upstrm.size());
                for (size_t i = 0; i < upstrm.size(); i++)
                        it[i] = upstrm[i].begin();

                for (size_t j = 0; j < path.size(); j++) {
                        BiNodeMap<int> nodeCounts;
                        for (size_t i = 0; i < upstrm.size(); i++) {
                                if (it[i] == upstrm[i].end()) {
                                        continue;
                                } else if (it[i]->first != path[j]) {
                                        it[i] = upstrm[i].end();
                                } else
                                        nodeCounts[next[i]] += (it[i]++)->second;
                        }

                        // all paths have reached the end
                        if (nodeCounts.empty())
                                return 0;

                        // compute conflict degree
                        int sumVal = 0, maxVal = -1; NodeID maxKey = 0;
                        for (const auto [id, count] : nodeCounts) {
                                sumVal += count;
                                if (count > maxVal) {
                                        maxVal = count;
                                        maxKey = id;
                                }
                        }

                        // check maximum conflict degree
                        if (CD(maxVal, sumVal - maxVal) <= MCD)
                                return maxKey;
                }

                return 0;
        }

        /**
         * Print all raw information
         */
        void print() const {
                for (size_t i = 0; i < upstrm.size(); i++)
                        std::cout << upstrm[i].getReverse() << " ---> " << next[i] << "\n";
        }
};


#endif
