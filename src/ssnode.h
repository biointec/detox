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

#ifndef SSNODE_H
#define SSNODE_H

#include "kmer/tkmer.h"
#include "dsnode.h"

// ============================================================================
// NODE REPRESENTATIVE CLASS
// ============================================================================

/**
 * Represents a node in a graph as a nodeID. Takes into account that
 * identifiers nodeID/-nodeID represent the same physical node
 */
class NodeRep {

private:
        NodeID nodeID;          // node identifier

public:
        /**
         * Default constructor
         */
        NodeRep() : nodeID(0) {}

        /**
         * Default constructor
         * @param target node identifier
         */
        NodeRep(NodeID target) : nodeID(abs(target)) {}

        /**
         * Get the node identifier
         * @return the node identifier
         */
        NodeID getNodeID() const {
                return nodeID;
        }

        /**
         * Operator == overloading
         * @return true or false
         */
        bool operator==(const NodeRep& rhs) const {
                return (nodeID == rhs.nodeID);
        }

        /**
         * Operator < overloading
         * @return true or false
         */
        bool operator<(const NodeRep& rhs) const {
                return nodeID < rhs.nodeID;
        }

        /**
         * Implicit conversion back to NodeID
         * @return Node identifier
         */
        operator NodeID() const {
                return nodeID;
        }

        /**
         * Compute a hash function
         * @return Hash value
         */
        size_t getHash() const
        {
                size_t w = nodeID;
                size_t hash = 0;
                w = ~w + (w << 21);             // key = (key << 21) - key - 1;
                w = w ^ (w >> 24);
                w = (w + (w << 3)) + (w << 8);  // key * 265
                w = w ^ (w >> 14);
                w = (w + (w << 2)) + (w << 4);  // key * 21
                w = w ^ (w >> 28);
                w = w + (w << 31);
                hash = hash ^ size_t(w);

                return hash;
        }
};

// ============================================================================
// HASH FUNCTION
// ============================================================================

struct NodeHash {
        size_t operator()(NodeRep nr) const {
                return nr.getHash();
        }
};

struct BiNodeHash {
        size_t operator()(NodeID id) const {
                size_t w = id;
                size_t hash = 0;
                w = ~w + (w << 21);             // key = (key << 21) - key - 1;
                w = w ^ (w >> 24);
                w = (w + (w << 3)) + (w << 8);  // key * 265
                w = w ^ (w >> 14);
                w = (w + (w << 2)) + (w << 4);  // key * 21
                w = w ^ (w >> 28);
                w = w + (w << 31);
                hash = hash ^ size_t(w);

                return hash;
        }
};

// ============================================================================
// NODE MAP
// ============================================================================

template<class T>
using NodeMap = google::sparse_hash_map<NodeRep, T, NodeHash>;

template<class T>
using BiNodeMap = google::sparse_hash_map<NodeID, T, BiNodeHash>;

template <typename T>
typename T::mapped_type get(T const& map, typename T::key_type const& key)
{
    typename T::const_iterator iter(map.find(key));
    return iter != map.end() ? iter->second : typename T::mapped_type();
}

// ============================================================================
// SINGLE STRANDED NODE CLASS
// ============================================================================

class SSNode {

private:
        DSNode *dsNode;         // reference to the double stranded node
        NodeID nodeID;          // identifier of the node

public:
        /**
         * Default constructor
         */
        SSNode() : nodeID(0), dsNode(NULL) {}

        /**
         * Constructor
         * @param dsNode Double stranded node
         * @param nodeID Identifier of the node
         */
        SSNode(DSNode* dsNode, NodeID nodeID) : dsNode(dsNode), nodeID(nodeID) {
                assert(nodeID != 0);
        }

        /**
         * Set the flag
         * @param flag True of false
         */
        void setFlag1(bool flag) {
                dsNode->setFlag1(flag);
        }

        /**
         * Get the flag
         * @return True of false
         */
        bool getFlag1() const {
                return dsNode->getFlag1();
        }

        /**
         * Set the flag
         * @param flag True of false
         */
        void setFlag2(bool flag) {
                dsNode->setFlag2(flag);
        }

        /**
         * Get the flag
         * @return True of false
         */
        bool getFlag2() const {
                return dsNode->getFlag2();
        }

        /**
         * Set the coverage
         * @param target The coverage
         */
        void setCov(Coverage target) {
                dsNode->setCov(target);
        }

        /**
         * Get the coverage
         * @return The coverage
         */
        Coverage getCov() const {
                return dsNode->getCov();
        }

        /**
         * Atomically increment the coverage
         * @param rhs Right hand side (default = 1.0f)
         */
        void incCov(Coverage rhs = 1.0f) {
                dsNode->incCov(rhs);
        }

        /**
         * Get the avarge kmer coverage
         * @return The average kmer coverage
         */
        Coverage getAvgCov() const {
                size_t ml = dsNode->getMarginalLength();
                return (ml > 0) ? dsNode->getCov() / (Coverage)ml : Coverage(0);
        }

        /**
         * Invalidate this node
         */
        void invalidate() {
                dsNode->invalidate();
        }

        /**
         * Check if a node is invalidated
         * @return True or false
         */
        bool isValid() const {
                if (nodeID == 0)
                        return false;
                return dsNode->isValid();
        }

        /**
         * Get the identifier of this node
         */
        NodeID getNodeID() const {
                return nodeID;
        }

        /**
         * Operator '==' overloading
         * @param rhs Right hand side SSNode
         * @return True if they're equal
         */
        bool operator==(const SSNode &rhs) const {
                if (dsNode != rhs.dsNode)
                        return false;
                return (nodeID == rhs.nodeID);
        }

        /**
         * Operator '!=' overloading
         * @param rhs Right hand side kmer
         * @return True if they're different
         */
        bool operator!=(const SSNode &rhs) const {
                return !(*this == rhs);
        }

        /**
         * Get the length of the node (# nucleotides in DNA string)
         * @return The length of the node
         */
        NodeLength length() const {
                return dsNode->getLength();
        }

        /**
         * The marginal length == length - k + 1
         * @return The marginal length of the node
         */
        NodeLength getMarginalLength() const {
                return dsNode->getMarginalLength();
        }

        /**
         * Get the number of left arcs
         * @return The number of left arcs
         */
        int numLeftArcs() const {
                return (nodeID > 0) ?
                        dsNode->numLeftArcs() : dsNode->numRightArcs();
        }

        /**
         * Get the number of right arcs
         * @return The number of right arcs
         */
        int numRightArcs() const {
                return (nodeID > 0) ?
                        dsNode->numRightArcs() : dsNode->numLeftArcs();
        }

        /**
         * Get the number of left arcs
         * @param The number of left arcs
         */
        void setNumLeftArcs(int numArcs) const {
                if (nodeID > 0)
                        dsNode->setNumLeftArcs(numArcs);
                else
                        dsNode->setNumRightArcs(numArcs);
        }

        /**
         * Get the number of right arcs
         * @param The number of right arcs
         */
        void setNumRightArcs(int numArcs) const {
                if (nodeID > 0)
                        dsNode->setNumRightArcs(numArcs);
                else
                        dsNode->setNumLeftArcs(numArcs);
        }

        /**
         * Get an iterator pointing to the first left arc
         * @return an iterator pointing to the first left arc
         */
        ArcIt leftBegin() const {
                return (nodeID > 0) ?
                        dsNode->leftBegin(false) : dsNode->rightBegin(true);
        }

        /**
         * Get an iterator pointing past the last left arc
         * @retur nan iterator pointing to the last left arc
         */
        ArcIt leftEnd() const {
                return (nodeID > 0) ?
                        dsNode->leftEnd(false) : dsNode->rightEnd(true);
        }

        /**
         * Get an iterator pointing to the first right arc
         * @return an iterator pointing to the first left arc
         */
        ArcIt rightBegin() const {
                return (nodeID > 0) ?
                        dsNode->rightBegin(false) : dsNode->leftBegin(true);
        }

        /**
         * Get an iterator pointing past the last right arc
         * @retur nan iterator pointing to the last right arc
         */
        ArcIt rightEnd() const {
                return (nodeID > 0) ?
                        dsNode->rightEnd(false) : dsNode->leftEnd(true);
        }

        /**
         * Delete the left arcs
         */
        void deleteAllLeftArcs() {
                if (nodeID > 0)
                        dsNode->deleteLeftArcs();
                else
                        dsNode->deleteRightArcs();
        }

        /**
         * Delete the right arcs
         */
        void deleteAllRightArcs() {
                if (nodeID > 0)
                        dsNode->deleteRightArcs();
                else
                        dsNode->deleteLeftArcs();
        }

        /**
         * Get the pointer to the left right arc
         * @return The pointer to the left right arc
         */
        Arc* getFirstLeftArc() {
                if (nodeID > 0)
                        return dsNode->getFirstLeftArc();
                return dsNode->getFirstRightArc();
        }

        /**
         * Get the pointer to the first right arc
         * @return The pointer to the first right arc
         */
        Arc* getFirstRightArc() {
                if (nodeID > 0)
                        return dsNode->getFirstRightArc();
                return dsNode->getFirstLeftArc();
        }

        /**
         * Set the pointer to the left right arc
         * @param target The pointer to the left right arc
         */
        void setFirstLeftArc(Arc* target) {
                if (nodeID > 0)
                        return dsNode->setFirstLeftArc(target);
                return dsNode->setFirstRightArc(target);
        }

        /**
         * Set the pointer to the first right arc
         * @param target The pointer to the first right arc
         */
        void setFirstRightArc(Arc* target) {
                if (nodeID > 0)
                        return dsNode->setFirstRightArc(target);
                return dsNode->setFirstLeftArc(target);
        }

        /**
         * Delete a specific left arc
         * @param targetID Identifier for the target node
         * @return True of the arc was deleted, false if the arc was not found
         */
        bool deleteLeftArc(NodeID targetID) {
                if (nodeID > 0)
                        return dsNode->deleteLeftArc(targetID);
                return dsNode->deleteRightArc(-targetID);
        }

        /**
         * Delete a specific right arc
         * @param targetID Identifier for the target node
         * @return True of the arc was deleted, false if the arc was not found
         */
        bool deleteRightArc(NodeID targetID) {
                if (nodeID > 0)
                        return dsNode->deleteRightArc(targetID);
                return dsNode->deleteLeftArc(-targetID);
        }

        /**
         * Get a specific left arc
         * @param targetID Identifier for the target node
         * @return Pointer to the specific arc, NULL if not found
         */
        Arc* leftArc(NodeID targetID) const {
                if (nodeID > 0)
                        return dsNode->leftArc(targetID);
                return dsNode->rightArc(-targetID);
        }

        /**
         * Get a specific right arc
         * @param targetID Identifier for the target node
         * @return Pointer to the specific arc, NULL if not found
         */
        Arc* rightArc(NodeID targetID) const {
                if (nodeID > 0)
                        return dsNode->rightArc(targetID);
                return dsNode->leftArc(-targetID);
        }

        /**
         * Get the nodeID a specific left arc
         * @param nucleotide Nucleotide under consideration
         * @return Node identifier
         */
       NodeID getLeftArcNodeID(char nucleotide) const {
                for (ArcIt it = leftBegin(); it != leftEnd(); it++) {
                        SSNode n(dsNode - abs(nodeID) + abs(it->getNodeID()),
                                 it->getNodeID());
                        if (n.peekNucleotideMarginalRight() == nucleotide)
                                return it->getNodeID();
                }

                return 0;
        }

        /**
         * Get a specific right arc
         * @param nucleotide Nucleotide under consideration
         * @return Pointer to the specific arc, NULL if not found
         */
        NodeID getRightArcNodeID(char nucleotide) const {
                for (ArcIt it = rightBegin(); it != rightEnd(); it++) {
                        SSNode n(dsNode - abs(nodeID) + abs(it->getNodeID()),
                                 it->getNodeID());
                        if (n.peekNucleotideMarginalLeft() == nucleotide)
                                return it->getNodeID();
                }

                return 0;
        }

        /**
         * Replace a left arc with a different one
         * @param origID Identifier for the original target node
         * @param newID Identifier for the new target node
         */
        void replaceLeftArc(NodeID origID, NodeID newID) {
                if (nodeID > 0) {
                        dsNode->leftArc(origID)->setNodeID(newID);
                } else {
                        dsNode->rightArc(-origID)->setNodeID(-newID);
                }
        }

        /**
         * Replace a right arc with a different one
         * @param origID Identifier for the original target node
         * @param newID Identifier for the new target node
         */
        void replaceRightArc(NodeID origID, NodeID newID) {
                if (nodeID > 0)
                        dsNode->rightArc(origID)->setNodeID(newID);
                else
                        dsNode->leftArc(-origID)->setNodeID(-newID);
        }

        /**
         * Get the sequence of this node
         * @return stl string containing the sequence
         */
        std::string getSequence() const {
                std::string seq = dsNode->getSequence();
                if (nodeID < 0)
                        Nucleotide::revCompl(seq);

                return seq;
        }

        /**
         * Get a subsequence of this node
         * @param pos Position of the first character to include
         * @param len Length of the subsequence
         * @return stl string containing the sequence
         */
        std::string substr(size_t pos = 0,
                           size_t len = std::string::npos) const {
                if (pos >= length())
                        return std::string();   // empty string
                len = std::min<size_t>(len, length() - pos);

                if (nodeID > 0)
                        return dsNode->substr(pos, len);
                else {
                        pos = length() - len - pos;
                        return Nucleotide::getRevCompl(dsNode->substr(pos, len));
                }
        }

        /**
         * Get a nucleotide at a specified position ('-' for out-of bounds)
         * @param pos Position in the sequence
         * @return Nucleotide at specified position
         */
        char getNucleotide(NodePosition pos) const {
                // check for out-of-bounds
                if (pos >= length())
                        return '-';
                if (nodeID < 0)
                        return Nucleotide::getComplement(dsNode->getNucleotide(length() - pos - 1));
                else
                        return dsNode->getNucleotide(pos);
        }

        /**
         * Set the sequence of this node
         * @param str String containing only 'A', 'C', 'G' and 'T'
         */
        void setSequence(const std::string& str) {
                if (nodeID > 0)
                        dsNode->setSequence(str);
                else
                        dsNode->setSequence(Nucleotide::getRevCompl(str));
        }

        /**
         * Get the left kmer of this node
         * @return The left kmer of this node
         */
        Kmer getLeftKmer() const {
                const std::string& seq = getSequence();
                return Kmer(seq);
        }

        /**
         * Get the right kmer of this node
         * @return The right kmer of this node
         */
        Kmer getRightKmer() const {
                const std::string& seq = getSequence();
                return Kmer(seq, seq.size() - Kmer::getK());
        }

        /**
         * Get the leftmost nucleotide of this node
         * @return The leftmost nucleotide
         */
        char peekNucleotideLeft() const {
                return (nodeID > 0) ? dsNode->peekNucleotideLeft() :
                        Nucleotide::getComplement(dsNode->peekNucleotideRight());
        }

        /**
         * Get the rightmost nucleotide of this node
         * @return The rightmost nucleotide
         */
        char peekNucleotideRight() const {
                return (nodeID > 0) ? dsNode->peekNucleotideRight() :
                        Nucleotide::getComplement(dsNode->peekNucleotideLeft());
        }

        /**
         * Get the nucleotide at position k - 1
         * @return The nucleotide at position k - 1
         */
        char peekNucleotideMarginalLeft() const {
                return (nodeID > 0) ? dsNode->peekNucleotideMarginalLeft() :
                        Nucleotide::getComplement(dsNode->peekNucleotideMarginalRight());
        }

        /**
         * Get the nucleotide at position size - k
         * @return The nucleotide at position size - k
         */
        char peekNucleotideMarginalRight() const {
                return (nodeID > 0) ? dsNode->peekNucleotideMarginalRight() :
                        Nucleotide::getComplement(dsNode->peekNucleotideMarginalLeft());
        }

        /**
         * Move the right arcs from another source node to the present node
         * @param src Source node
         */
        void inheritRightArcs(SSNode& src) {
                // make sure the node has currently no right arcs
                assert(numRightArcs() == 0);

                // update the arc information for the connected nodes
                for (ArcIt it = src.rightBegin(); it != src.rightEnd(); it++) {
                        SSNode r(dsNode - abs(nodeID) + abs(it->getNodeID()),
                                 it->getNodeID());
                        r.replaceLeftArc(src.getNodeID(), nodeID);
                }

                // move the right arcs
                setNumRightArcs(src.numRightArcs());
                src.setNumRightArcs(0);
                setFirstRightArc(src.getFirstRightArc());
                // if src and present node have opposite sign
                if ((nodeID < 0) != (src.getNodeID() < 0)) {
                        if (nodeID > 0)
                                dsNode->swapSignRightArcs();
                        else
                                dsNode->swapSignLeftArcs();
                }
        }
};

#endif
