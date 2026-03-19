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

#ifndef NODECHAIN_H
#define NODECHAIN_H

#include <vector>
#include <deque>
#include <set>

#include "ssnode.h"

// ============================================================================
// NODE CHAIN CLASS
// ============================================================================

class NodeChain : public std::vector<NodeID>
{
private:
        size_t count;   // number of times the nodechain was observed

public:
        /**
         * Default constructor
         */
        NodeChain() : count(0) {}

        /**
         * Constructor from an input vector
         * @param input Input vector
         * @param count Count
         */
        NodeChain(const std::vector<NodeID>& input, size_t count = 1) :
                count(count)
        {
                for (const auto& it : input)
                        push_back(it);
        }

        /**
         * Constructor from an input vector range [first, last[
         * @param first Iterator the first element
         * @param last Iterator past the last element
         * @param count Count
         */
        NodeChain(std::vector<NodeID>::const_iterator first,
                  std::vector<NodeID>::const_iterator last, size_t count = 1) :
                        count(count)
        {
                if (last <= first)
                        return;
                reserve(last - first);
                for (auto it = first; it < last; it++)
                        push_back(*it);
        }

        /**
         * Increment the count by val
         * @param val Value (default = 1)
         */
        void incrCount(size_t val = 1) {
                count += val;
        }

        /**
         * Get the count
         * @return count
         */
        size_t getCount() const {
                return count;
        }

        /**
         * Set the count
         * @param target Target count
         */
        void setCount(size_t target) {
                count = target;
        }

        /**
         * Operator == overloading (disregard count value)
         * @param rhs Right hand side
         */
        bool operator==(const NodeChain& rhs) const {
                if (size() != rhs.size())
                        return false;
                for (size_t i = 0; i < size(); i++)
                        if ((*this)[i] != rhs[i])
                                return false;
                return true;
        }

        /**
         * Reverse the node chain
         */
        void reverse() {
                std::reverse(begin(), end());
        }

        /**
         * Get the reverse of the node chain
         * @return The reverse of the node chain
         */
        NodeChain getReverse() const {
                NodeChain copy = *this;
                copy.reverse();
                return copy;
        }

        /**
         * Reverse-complement the node chain
         */
        void revCompl() {
                reverse();
                for (auto& e : *this)
                        e = -e;
        }

        /**
         * Get the reverse complement of the node chain
         * @return The reverse complementary chain
         */
        NodeChain getRevCompl() const {
                NodeChain copy = *this;
                copy.revCompl();
                return copy;
        }

        /**
         * Get the representative node chain
         * @return The representative node chain
         */
        NodeChain getRepresentative() const {
                NodeChain RC = getRevCompl();
                return (RC < *this) ? RC : *this;
        }

        bool contains(const std::set<NodeRep>& noi) const {
                for (auto it : *this)
                        if (noi.find(NodeRep(it)) != noi.end())
                                return true;
                return false;
        }

        bool contains(NodeRep nr) const {
                for (auto it : *this)
                        if (NodeRep(it) == nr)
                                return true;
                return false;
        }

        bool contains(NodeID nodeID) const {
                for (auto it : *this)
                        if (it == nodeID)
                                return true;
                return false;
        }

        bool isRepeated(NodeRep nr) const {
                int count = 0;
                for (auto it : *this)
                        if (NodeRep(it) == nr)
                                count++;
                return (count > 1);
        }

        bool isRepeated(NodeID nodeID) const {
                int count = 0;
                for (auto it : *this)
                        if (it == nodeID)
                                count++;
                return (count > 1);
        }

        /**
         * Operator<< overloading
         * @param out Output file stream (input)
         * @param ncc NodeChainContainer to display
         * @return Output file stream
         */
        friend std::ostream &operator<<(std::ostream &out, const NodeChain& nc);

        static bool compCov(const NodeChain& l, const NodeChain& r) {
                return (l.getCount() < r.getCount());
        }

        static bool covHiToLo(const NodeChain& l, const NodeChain& r) {
                return (l.getCount() > r.getCount());
        }

        /**
         * @brief Sort short to long, then high count to low count
         * @param l Left hand side
         * @param r Right hand side
         * @return true or false
         */
        static bool compH2LLex(const NodeChain& l, const NodeChain& r) {
                if (l.getCount() != r.getCount())
                        return (l.getCount() > r.getCount());
                for (size_t i = 0; (i < l.size()) && (i < r.size()); i++)
                        if (l[i] != r[i])
                                l[i] < r[i];
                return l.size() < r.size();
        }

        /**
         * @brief Sort short to long, then high count to low count
         * @param l Left hand side
         * @param r Right hand side
         * @return true or false
         */
        static bool compS2LH2L(const NodeChain& l, const NodeChain& r) {
                if (l.size() != r.size())
                        return (l.size() < r.size());
                return (l.getCount() > r.getCount());
        }

        /**
         * @brief Sort long to short, then high count to low count
         * @param l Left hand side
         * @param r Right hand side
         * @return true or false
         */
        static bool compL2SH2L(const NodeChain& l, const NodeChain& r) {
                if (l.size() != r.size())
                        return (l.size() > r.size());
                return (l.getCount() > r.getCount());
        }

        /**
         * Sort by length (long to short), then lexigraphically
         */
        static bool compLenLex(const NodeChain& l, const NodeChain& r) {
                if (l.size() != r.size())
                        return (l.size() > r.size());
                for (size_t i = 0; i < l.size(); i++)
                        if (l[i] != r[i])
                                l[i] < r[i];
                return false;
        }

        /**
         * Merge two overlapping node chains
         * @param L Left node chain
         * @param R Right node chain
         * @return Merged node chain, empty if L and R cannot be merged
         */
        static NodeChain merge(const NodeChain& L, const NodeChain& R) {
                if (L.empty())
                        return R;
                if (R.empty())
                        return L;

                for (size_t s = 0; s < L.size(); s++) {
                        bool canMerge = true;
                        for (size_t i = s, j = 0;
                             i < L.size() && j < R.size(); i++, j++) {
                                if (L[i] == R[j])
                                        continue;
                                canMerge = false;
                                break;
                        }

                        if (!canMerge)
                                continue;

                        NodeChain m;
                        for (size_t i = 0; i < s; i++)
                                m.push_back(L[i]);
                        for (size_t j = 0; j < R.size(); j++)
                                m.push_back(R[j]);
                        m.setCount(std::min(L.getCount(), R.getCount()));
                        return m;
                }
                return NodeChain();     // return empty node chain
        }
};

// ============================================================================
// CONSENSUS NODE CHAIN CLASS
// ============================================================================

class ConsNodeChain : public std::deque<std::pair<NodeID, int>>
{
public:
        // inherit base constructors (C++11)
        using std::deque<std::pair<NodeID, int>>::deque;

        /**
         * Construct a Consensus Node Chain from a regulator Node Chain
         * @param nc Regular node chain
         */
        ConsNodeChain(const NodeChain& nc) {
                int c = nc.getCount();
                std::transform(nc.begin(), nc.end(), std::back_inserter(*this),
                        [c](const NodeID& id){ return std::make_pair(id, c); });
        }

        /**
         * Operator<< overloading
         * @param out Output file stream (input)
         * @param ncc NodeChainContainer to display
         * @return Output file stream
         */
        friend std::ostream &operator<<(std::ostream &out,
                                        const ConsNodeChain& nc);

        /**
         * Reverse the consensus node chain
         */
        void reverse() {
                std::reverse(begin(), end());
        }

        /**
         * Get the reverse of the node chain
         * @return The reverse of the node chain
         */
        ConsNodeChain getReverse() const {
                ConsNodeChain copy = *this;
                copy.reverse();
                return copy;
        }

        /**
         * Reverse complement the consensus node chain
         */
        void revCompl() {
                reverse();
                for (auto& [n, w] : *this)
                        n = -n;
        }

        /**
         * Get the reverse complement of the consensus node chain
         * @return The reverse complement consensus node chain
         */
        ConsNodeChain getRevCompl() const {
                ConsNodeChain copy = *this;
                copy.revCompl();
                return copy;
        }

        /**
         * Check if the consensus node chain is prefixed by a node chain
         * @param nc Node chain to act as a prefix
         * @return True if mergeable, false otherwise
         */
        bool isPrefixed(const NodeChain& nc) const {
                if (nc.size() > size())
                        return false;           // node chain too long
                for (size_t i = 0; i < nc.size(); i++)
                        if ((*this)[i].first != nc[i])
                                return false;   // not mergeable
                return true;
        }

        /**
         * Merge a node chain into a consensus node chain
         * @param nc Node chain to merge
         * @return True if mergeable, false otherwise
         */
        bool mergeNodeChain(const NodeChain& nc) {
                // check if the node chain is mergeable
                for (size_t i = 0; i < std::min<size_t>(size(), nc.size()); i++)
                        if ((*this)[i].first != nc[i])
                                return false;   // not mergeable
                // do the actual merging
                const int c = nc.getCount();
                size_t minLen = std::min<size_t>(size(), nc.size());
                for (size_t i = 0; i < minLen; i++)
                        (*this)[i].second += c;
                std::transform(nc.begin() + minLen, nc.end(),
                               std::back_inserter(*this),
                               [c](const NodeID& id){ return std::make_pair(id, c); });
                return true;
        }

        /**
         * Merge a consensus node chain into a consensus node chain
         * @param nc Consensus node chain to merge
         * @return True if mergeable, false otherwise
         */
        bool mergeNodeChain(const ConsNodeChain& nc) {
                // check if the consensus node chain is mergeable
                for (size_t i = 0; i < std::min<size_t>(size(), nc.size()); i++)
                        if ((*this)[i].first != nc[i].first)
                                return false;   // not mergeable
                // do the actual merging
                for (size_t i = 0; i < std::min<size_t>(size(), nc.size()); i++)
                        (*this)[i].second += nc[i].second;

                if (size() < nc.size())
                        insert(end(), nc.begin() + size(), nc.end());
                return true;
        }

        /**
         * Convert a consensus node chain to a regular node chain
         * @return Regular node chain
         */
        NodeChain toNodeChain() const
        {
                if (empty())
                        return NodeChain();

                std::vector<NodeID> path(size());
                int weight = front().second;
                for (size_t i = 0; i < size(); i++) {
                        path[i] = (*this)[i].first;
                        weight = std::min<int>((*this)[i].second, weight);
                }

                return NodeChain(path.begin(), path.end(), weight);
        }

        /**
         * Does the consensus node chain contain a certain node?
         * @param nodeID Node indentifier
         * @return true or false
         */
        bool contains(NodeID nodeID) const {
                for (const auto& it : *this)
                        if (it.first == nodeID)
                                return true;
                return false;
        }

        /**
         * Does the consensus node chain contain a certain node?
         * @param nodeID Node indentifier
         * @return true or false
         */
        bool containsRep(NodeRep nr) const {
                for (const auto& it : *this)
                        if (NodeRep(it.first) == nr)
                                return true;
                return false;
        }

        /**
         * Merge two overlapping consensus node chains. The resulting
         * consensus node chain is guaranteed to be a path in the graph.
         * @param X First consensus node chain
         * @param Y Second consensus node chain
         * @param MCD Maximum conflict degree
         * @return Merged node chain, empty if X and Y cannot be merged
         */
        static ConsNodeChain merge(const ConsNodeChain& L,
                                   const ConsNodeChain& R,
                                   float MCD);
};

#endif
