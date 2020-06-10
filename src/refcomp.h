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

#ifndef REFCOMP_H
#define REFCOMP_H

#include "global.h"
#include "kmernpp.h"

#include <string>
#include <vector>
#include <algorithm>

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class DBGraph;
class Settings;
class WorkLoadBalancer;
class GraphAligner;

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
         */
        NodeChain(const std::vector<NodeID>& input) : count(1) {
                for (const auto& it : input)
                        push_back(it);
        }

        /**
         * Increment the count by one
         */
        void incrCount() {
                count++;
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
         * Operator < overloading
         * @param rhs Right hand side
         */
        bool operator<(const NodeChain& rhs) const {
                if (size() != rhs.size())
                        return size() < rhs.size();
                for (size_t i = 0; i < size(); i++)
                        if ((*this)[i] != rhs[i])
                                return (*this)[i] < rhs[i];
                return false;
        }

        /**
         * Operator < overloading
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
         * Get the reverse complement of the node chain
         * @return The reverse complementary chain
         */
        NodeChain getReverseComplement() const {
                NodeChain copy = *this;

                std::reverse(copy.begin(), copy.end());
                for (size_t i = 0; i < copy.size(); i++)
                        copy[i] = -copy[i];

                return copy;
        }

        /**
         * Get the representative node chain
         * @return The representative node chain
         */
        NodeChain getRepresentative() const {
                NodeChain RC = getReverseComplement();
                return (RC < *this) ? RC : *this;
        }

        /**
         * Operator<< overloading
         * @param out Output file stream (input)
         * @param ncc NodeChainContainer to display
         * @return Output file stream
         */
        friend std::ostream &operator<<(std::ostream &out, const NodeChain& nc);
};

// ============================================================================
// ALIGNED SEGMENT CLASS
// ============================================================================

/*
 * The alignment of a sequence (e.g. reference genome) to a de Bruijn graph is
 * represented as a number of consecutive alignment segments:
 * a) Contigs: subsequences that align perfectly to the dBG
 * b) Parallel: subsequences of the sequence that align to a parallel path
 * c) Break: subsequences that do not align to the dBG
 * d) Deletion: subsequences that are missing in the dBG
 * e) Insertion: subsequences that only exist in the dBG
 * All positions are represented in k-mer space
 */

enum class AlnSegmentType { CONTIG, PARALLEL, BREAK, DELETION, INSERTION };

class AlnSegment {

public:
        AlnSegmentType type;    // type of segment

        size_t seqBegin;        // begin position in the subsequence
        NodeID srcNodeID;       // node ID of the first node
        size_t srcNodeBegin;    // begin position in the start node

        size_t seqEnd;          // end position in the subsequence
        NodeID dstNodeID;       // node ID of the last node
        size_t dstNodeEnd;      // end position in the last node

        /**
         * Default constructor
         * @param type Type of segment (defaults to CONTIG)
         */
        AlnSegment(AlnSegmentType type = AlnSegmentType::CONTIG) : type(type),
                seqBegin(0), srcNodeID(0), srcNodeBegin(0), seqEnd(0),
                dstNodeID(0), dstNodeEnd(0) {};

        /**
         * Set the start point of an aligned segment
         * @param seqBegin Begin position in the sequence
         * @param srcNodeID Node ID of the first node
         * @param srcNodeBegin Begin position in the start node
         */
        void setBegin(size_t seqBegin, NodeID srcNodeID, size_t srcNodeBegin) {
                this->seqBegin = seqBegin;
                this->srcNodeID = srcNodeID;
                this->srcNodeBegin = srcNodeBegin;
        }

        /**
         * Set the end point of an aligned segment
         * @param seqEnd End position in the sequence
         * @param dstNodeID Node ID of the last node
         * @param dstNodeEnd End position of the last node
         */
        void setEnd(size_t seqEnd, NodeID dstNodeID, size_t dstNodeEnd) {
                this->seqEnd = seqEnd;
                this->dstNodeID = dstNodeID;
                this->dstNodeEnd = dstNodeEnd;
        }

        /**
         * Get the length of segment
         * @return The length of the sequence
         */
        size_t getSeqLength() const {
                return seqEnd - seqBegin;
        }

        /**
         * Get a node-position pair of the begin position of the segment
         * @return node-position pair
         */
        NodePosPair getBegin() const {
                return NodePosPair(srcNodeID, srcNodeBegin);
        }

        /**
         * Get a node-position pair of the end position of the segment
         * @return node-position pair
         */
        NodePosPair getEnd() const {
                return NodePosPair(dstNodeID, dstNodeEnd);
        }

        /**
         * Operator << overloading
         * @param out Output stream
         * @param s Segment to display
         * @return Output stream
         */
        friend std::ostream &operator<<(std::ostream &out, const AlnSegment &s);
};

// ============================================================================
// REFERENCE COMPARISON CLASS
// ============================================================================

class RefComp {

private:
        const DBGraph& dBG;                     // const-reference to dBG
        const Settings& settings;               // const-reference to settings
        std::vector<std::string> sequence;      // reference sequences

        /**
         * Read a fasta file from disk and store the sequences
         * @param fastaFilename Filename of fasta file with sequences
         */
        void readFastaFile(const std::string& fastaFilename);

        /**
         * Align a sequence to de Bruijn graph
         * @param seq Sequence to align
         * @param segment Alignment segment vector (output)
         */
        void alignSequence(const std::string& seq,
                           std::vector<AlnSegment>& segment) const;

        /**
         * @brief Annotate a single breakpoint
         * @param seq Sequence to align
         * @param segment Alignment segment vector (output)
         * @param graphAln Graph aligner object
         */
        void annotateBreakpoint(size_t i, const std::string& seq,
                                std::vector<AlnSegment>& segment,
                                GraphAligner& graphAln) const;

        /**
         * @brief Annotate the breakpoints in the alignment (thread entry)
         * @param myID Thread identifier
         * @param wlb Workload balancer
         * @param seq Sequence to align
         * @param segment Alignment segment vector (output)
         */
        void annotateBreakpointThread(size_t myID, WorkLoadBalancer& wlb,
                                      const std::string& seq,
                                      std::vector<AlnSegment>& segment) const;

        /**
         * @brief Annotate the breakpoints in the alignment
         * @param seq Sequence to align
         * @param segment Alignment segment vector (output)
         * @param output Progress indicator output
         */
        void annotateBreakpoints(const std::string& seq,
                                 std::vector<AlnSegment>& segment,
                                 const std::string& output) const;

        /**
         * @brief Output segments to the screen
         * @param seq Sequence to align
         * @param segment Segments to output
         */
        void printSegments(const std::string& seq,
                           const std::vector<AlnSegment>& segment) const;

        /**
         * Find k-mer in the Kmernodetable
         * @param kmer k-mer respresented as a string (can contain non-ACTG)
         * @return < NodeID, NodePosition > pair
         */
        /*NodePosPair findNPP(const std::string& kmer) const {
                assert(kmer.size() == Kmer::getK());

                // if one of the characters is non-ACTG return "not found"
                for (char c : kmer)
                        if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
                                return NodePosPair(0, 0);

                return table.find(Kmer(kmer));
        }*/

public:
        /**
         * Constructor
         * @param dBG Const-ref to de Bruijn graph
         */
        RefComp(const DBGraph& dBG, const Settings& settings, const std::string& fastaFilename) :
                dBG(dBG), settings(settings)
        {
                readFastaFile(fastaFilename);
        }

        /**
         * Align reference sequences to a de Bruijn graph
         */
        void alignSequences();

        /**
         * Get the true node chains from the reference sequence
         * @param nodeChain Node chains (output)
         */
        void getTrueNodeChain(std::vector<NodeChain>& nodeChain);

        /**
         * Calculate the true node multiplicity
         * @param table Populated <kmer, NodePosPair> table
         * @param nodeMult Node multiplicity vector (output)
         * @param edgeMult Edge multiplicity vector (output)
         */
        void getTrueMultiplicity(const KmerNPPTable& table,
                                 std::vector<int>& nodeMult,
                                 std::vector<int>& edgeMult);
};

#endif
