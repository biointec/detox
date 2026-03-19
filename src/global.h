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

#ifndef GLOBAL_H
#define GLOBAL_H

#include <cstdlib>
#include <cstdint>
#include <utility>

#define FLOAT_SMALL 1E-6
#define DOUBLE_SMALL 1E-15

// The number of recordBlocks simulataneously held in memory. Higher values
// result in more parallel chunks at the cost of increased memory use.
#define NUM_RECORD_BLOCKS 2

// Tip clipping
#define MAX_TIP_LEN 150
#define MAX_TIP_INDELS 3
#define MAX_TIP_ITERATIONS 100

// Bubble popping
#define MAX_BUBBLE_LEN 150
#define MAX_BUBBLE_ED_PER_KMER 5
#define MAX_BUBBLE_INDELS 3
#define MAX_BUBBLE_ITERATIONS 100

#define MAX_CHIM_LEN 150
#define MIN_TIP_MERGE_LEN 250

#define INDEP_OBS 250

// Annotate breakpoint
#define MAX_NODES_REWIND 3
#define MAX_NUCLEOTIDES_REWIND 100

// size of the ends of a contig used for paired-end read anchoring (scaffolding)
#define MAX_CONTIG_ANCHOR_LEN 150

// ============================================================================
// TYPEDEFS
// ============================================================================

typedef int32_t NodeID;                 // max 2 billion nodes
typedef int32_t ArcID;                  // max 2 billion arcs
typedef int32_t ContigID;
typedef uint64_t NucleotideID;
typedef uint32_t NodeLength;            // max length is 4 billion
typedef NodeLength NodePosition;        // position in a contig or read
typedef float Coverage;

template<size_t numBytes>
class TKmer;

template<size_t numBytes>
struct TKmerHash;

typedef std::pair<NodeID, NodeID> NodePair;

#endif
