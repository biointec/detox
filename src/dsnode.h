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

#ifndef DSNODE_H
#define DSNODE_H

#include "arc.h"
#include "kmer/tstring.h"

#include <map>
#include <set>
#include <atomic>

// ============================================================================
// DOUBLE STRANDED NODE CLASS
// ============================================================================

class DSNode {

private:
        typedef union {
                struct Packed {
                        uint8_t numLeft:3;              // [0...4]
                        uint8_t numRight:3;             // [0...4]
                        uint8_t invalid:1;
                        uint8_t flag:1;
                } p;
                uint8_t up;
        } Bitfield;

        TString sequence;       // DNA sequence
        Arc* lArc;              // Pointer to the first left arc
        Arc* rArc;              // Pointer to the first right arc
        Bitfield arcInfo;       // number of arcs at each node

        std::atomic<Coverage> cov;
        std::atomic<bool> myFlag1;
        std::atomic<bool> myFlag2;

public:
        /**
         * Default constructor
         */
        DSNode() : lArc(NULL), rArc(NULL), cov(0), myFlag1(false), myFlag2(false)
        {
                arcInfo.up = 0;
        }

        DSNode& operator=(DSNode&& rhs) {

                if (this == &rhs)
                        return *this;

                // copy/move all members
                sequence = std::move(rhs.sequence);
                lArc = rhs.lArc;
                rArc = rhs.rArc;
                arcInfo = rhs.arcInfo;
                cov = rhs.cov.load();
                myFlag1 = rhs.myFlag1.load();
                myFlag2 = rhs.myFlag2.load();

                // load rhs in an empty state
                rhs.lArc = rhs.rArc = NULL;
                rhs.arcInfo.up = 0;
                rhs.cov = 0;
                rhs.myFlag1 = false;
                rhs.myFlag2 = false;

                return *this;
        }

        /**
         * Get the pointer to the left right arc
         * @return The pointer to the left right arc
         */
        Arc* getFirstLeftArc() const {
                return lArc;
        }

        /**
         * Get the pointer to the first right arc
         * @return The pointer to the first right arc
         */
        Arc* getFirstRightArc() const {
                return rArc;
        }

        /**
         * Set the pointer to the left right arc
         * @param target The pointer to the left right arc
         */
        void setFirstLeftArc(Arc* target) {
                lArc = target;
        }

        /**
         * Set the pointer to the first right arc
         * @param target The pointer to the first right arc
         */
        void setFirstRightArc(Arc* target) {
                rArc = target;
        }

        /**
         * Set the coverage
         * @param target The coverage
         */
        void setCov(Coverage target) {
                cov = target;
        }

        /**
         * Get the coverage
         * @return The coverage
         */
        Coverage getCov() const {
                return cov;
        }

        /**
         * Atomically increment the coverage
         * @param rhs Right hand side (default = 1.0f)
         */
        void incCov(Coverage rhs = 1.0f) {
                auto current = cov.load();
                while (!cov.compare_exchange_weak(current, current + rhs));
        }

        /**
         * Set the flag
         * @param flag True of false
         */
        void setFlag1(bool flag) {
                myFlag1 = flag;
        }

        /**
         * Get the flag
         * @return True of false
         */
        bool getFlag1() const {
                return myFlag1;
        }

        /**
         * Set the flag
         * @param flag True of false
         */
        void setFlag2(bool flag) {
                myFlag2 = flag;
        }

        /**
         * Get the flag
         * @return True of false
         */
        bool getFlag2() const {
                return myFlag2;
        }

        /**
         * Swap the sign of the right arcs
         */
        void swapSignRightArcs() {
                for (int i = 0; i < arcInfo.p.numRight; i++)
                        rArc[i].setNodeID(-rArc[i].getNodeID());
        }

        /**
         * Swap the sign of the left arcs
         */
        void swapSignLeftArcs() {
                for (int i = 0; i < arcInfo.p.numLeft; i++)
                        lArc[i].setNodeID(-lArc[i].getNodeID());
        }

        /**
         * Invalidate a node (= mark as deleted)
         */
        void invalidate() {
                arcInfo.p.invalid = 1;
                sequence.clear();
        }

        /**
         * Check if a node is valid
         * @return True or false
         */
        bool isValid() const {
                return (arcInfo.p.invalid == 0);
        }

        /**
         * Get the length of the node (# nucleotides in DNA string)
         * @return The length of the node
         */
        NodeLength getLength() const {
                return sequence.getLength();
        }

        /**
         * The marginal length == length - k + 1
         * @return The marginal length of the node
         */
        NodeLength getMarginalLength() const {
                size_t l = getLength();
                return (l >= Kmer::getK()) ? l - Kmer::getK() + 1 : 0;
        }

        /**
         * Set the number of left arcs
         * @param numLeft The number of left arcs
         */
        void setNumLeftArcs(int numLeft) {
                arcInfo.p.numLeft = (uint8_t)numLeft;
        }

        /**
         * Set the number of right arcs
         * @param numright The number of right arcs
         */
        void setNumRightArcs(int numRight) {
                arcInfo.p.numRight = (uint8_t)numRight;
        }

        /**
         * Get the number of left arcs
         * @return The number of left arcs
         */
        int numLeftArcs() const {
                return arcInfo.p.numLeft;
        }

        /**
         * Get the number of right arcs
         * @return The number of right arcs
         */
        int numRightArcs() const {
                return arcInfo.p.numRight;
        }

        /**
         * Delete all left arcs
         */
        void deleteLeftArcs() {
                for (int i = 0; i < arcInfo.p.numLeft; i++)
                        lArc[i].deleteArc();
                arcInfo.p.numLeft = 0;
        }

        /**
         * Delete all right arcs
         */
        void deleteRightArcs() {
                for (int i = 0; i < arcInfo.p.numRight; i++)
                        rArc[i].deleteArc();
                arcInfo.p.numRight = 0;
        }

        /**
         * Get a specific left arc
         * @param nodeID Identifier for the target node
         * @return Pointer to the specific arc, NULL if not found
         */
        Arc* leftArc(NodeID nodeID) {
                for (int i = 0; i < arcInfo.p.numLeft; i++)
                        if (lArc[i].getNodeID() == nodeID)
                                return lArc + i;
                return NULL;
        }

        /**
         * Get a specific right arc
         * @param nodeID Identifier for the target node
         * @return Pointer to the specific arc, NULL if not found
         */
        Arc* rightArc(NodeID nodeID) {
                for (int i = 0; i < arcInfo.p.numRight; i++)
                        if (rArc[i].getNodeID() == nodeID)
                                return rArc + i;
                return NULL;
        }

        /**
         * Delete a specific left arc
         * @param targetID Identifier for the target node
         * @return True of the arc was deleted, false if the arc was not found
         */
        bool deleteLeftArc(NodeID targetID);

        /**
         * Delete a specific right arc
         * @param targetID Identifier for the target node
         * @return True of the arc was deleted, false if the arc was not found
         */
        bool deleteRightArc(NodeID targetID);

        /**
         * Get an iterator pointing to the first left arc
         * @param reversed True if you want arcs from the negative node
         * @return An iterator pointing to the first left arc
         */
        ArcIt leftBegin(bool reversed = false) const {
                return ArcIt(lArc, reversed);
        }

        /**
         * Get an iterator pointing past the last left arc
         * @param reversed True if you want arcs from the negative node
         * @return An iterator pointing to the last left arc
         */
        ArcIt leftEnd(bool reversed = false) const {
                return ArcIt(lArc + arcInfo.p.numLeft, reversed);
        }

        /**
         * Get an iterator pointing to the first right arc
         * @param reversed True if you want arcs from the negative node
         * @return An iterator pointing to the first left arc
         */
        ArcIt rightBegin(bool reversed = false) const {
                return ArcIt(rArc, reversed);
        }

        /**
         * Get an iterator pointing past the last right arc
         * @param reversed True if you want arcs from the negative node
         * @return An iterator pointing to the last right arc
         */
        ArcIt rightEnd(bool reversed = false) const {
                return ArcIt(rArc + arcInfo.p.numRight, reversed);
        }

        /**
         * Set the sequence of this node
         * @param str String containing only 'A', 'C', 'G' and 'T'
         */
        void setSequence(const std::string& str) {
                sequence.setSequence(str);
        }

        /**
         * Get the sequence of this node
         * @return The sequence of this node
         */
        std::string getSequence() const {
                return sequence.getSequence();
        }

        /**
         * Get a subsequence of this node
         * @param pos Position of the first character to include
         * @param len Length of the subsequence
         * @return stl string containing the sequence
         */
        std::string substr(size_t pos = 0,
                           size_t len = std::string::npos) const {
                return sequence.substr(pos, len);
        }

        /**
         * Get a nucleotide at a specified position ('-' for out-of bounds)
         * @param pos Position in the sequence
         * @return Nucleotide at specified position
         */
        char getNucleotide(NodePosition pos) const {
                // check for out-of-bounds
                if (pos >= getLength())
                        return '-';
                return sequence[pos];
        }

        /**
         * Get the tight sequence of this node
         * @return The tight sequence
         */
        const TString& getTSequence() const {
                return sequence;
        }

        /**
         * Get the leftmost nucleotide of this node
         * @return The leftmost nucleotide
         */
        char peekNucleotideLeft() const {
                return sequence.peekNucleotideLeft();
        }

        /**
         * Get the rightmost nucleotide of this node
         * @return The rightmost nucleotide
         */
        char peekNucleotideRight() const {
                return sequence.peekNucleotideRight();
        }

        /**
         * Get the nucleotide at position k - 1
         * @return The nucleotide at position k - 1
         */
        char peekNucleotideMarginalLeft() const {
                return sequence.peekNucleotideMarginalLeft();
        }

        /**
         * Get the nucleotide at position size - k
         * @return The nucleotide at position size - k
         */
        char peekNucleotideMarginalRight() const {
                return sequence.peekNucleotideMarginalRight();
        }

        /**
         * Get the leftmost kmer of this node
         * @return The leftmost kmer
         */
        Kmer getLeftKmer() const {
                return Kmer(sequence, 0);
        }

        /**
         * Get the rightmost kmer of this node
         * @return The rightmost kmer
         */
        Kmer getRightKmer() const {
                return Kmer(sequence, sequence.getLength() - Kmer::getK());
        }

        /**
         * Write a node to file
         * @param ofs Open output file stream
         */
        void write(std::ofstream& ofs, Arc* arcBase) const {
                ofs.write((char*)&cov,sizeof(cov));
                ArcID leftID = lArc - arcBase;
                ofs.write((char*)&leftID, sizeof(leftID));
                ArcID rightID = rArc - arcBase;
                ofs.write((char*)&rightID, sizeof(rightID));
                ofs.write((char*)&arcInfo, sizeof(arcInfo));
                sequence.write(ofs);
        }

        /**
         * Load a node from file
         * @param ifs Open input file stream
         */
        void read(std::ifstream& ifs, Arc* arcBase) {
                ifs.read((char*)&cov, sizeof(cov));
                ArcID leftID, rightID;
                ifs.read((char*)&leftID, sizeof(leftID));
                lArc = arcBase + leftID;
                ifs.read((char*)&rightID, sizeof(rightID));
                rArc = arcBase + rightID;
                ifs.read((char*)&arcInfo, sizeof(arcInfo));
                sequence.read(ifs);
        }
};

// ============================================================================
// HIGHER-ORDER NODE IDENTIFIER
// ============================================================================

typedef std::vector<NodeID> HoNodeID;
std::ostream& operator<<(std::ostream& out, const HoNodeID &id);

// ============================================================================
// HIGHER-ORDER EDGE IDENTIFIER
// ============================================================================

typedef std::pair<HoNodeID, HoNodeID> HoEdgeID;
std::ostream& operator<<(std::ostream& out, const HoEdgeID &id);

#endif
