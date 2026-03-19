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

#ifndef LIBRARY_H
#define LIBRARY_H

#include "util.h"
#include "ssnode.h"

#include <vector>
#include <string>
#include <filesystem>

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class DBGraph;
class NodeTracker;
enum ReadOrientation { UNKNOWN, FF, RF, FR };

// ============================================================================
// LIBRARY CLASS
// ============================================================================

class Library
{
private:
        typedef std::pair<std::string, std::string> StringPair;

        StringPair filename;            // if single-end then 2nd is empty
        StringPair baseFilename;        // if single-end then 2nd is empty
        float fragmentSize;             // fragment size
        float fragmentStd;              // fragment standard deviation

        int alnVersion;                 // alignment version

        /**
         * Infer the read orientation (FF, RF or FR)
         * @param dBG Const-reference to de Bruijn graph
         * @param numReads Number of reads to analyze
         * @return Paired-end read orientation
         */
        ReadOrientation inferReadOrientation(const DBGraph& dBG,
                                             size_t numReads = 100) const;

public:
        /**
         * Construct a new Library object
         * @param filename Input filename (if single-end then 2nd is empty)
         * @param baseFilename Base filename (if single-end then 2nd is empty)
         * @param fragmentSize Fragment size
         * @param fragmentStd Fragment standard deviation
         */
        Library(const StringPair& filename, const StringPair& baseFilename,
                float fragmentSize = 0.0f, float fragmentStd = 0.0f) :
                filename(filename), baseFilename(baseFilename),
                fragmentSize(fragmentSize), fragmentStd(fragmentStd),
                alnVersion(0) {};

        /**
         * Get the filename(s)
         * @return The filename(s) as a pair-of-strings
         */
        const StringPair& getFilename() const {
                return filename;
        }

        /**
         * Get the alignment filename(s)
         * @return The alignment filename(s) as a pair-of-strings
         */
        StringPair getAlnFilename() const {
                // strip path from read files
                namespace fs = std::filesystem;
                std::string base_fn_1 =
                        fs::path(baseFilename.first).filename().string();
                std::string base_fn_2 =
                        fs::path(baseFilename.second).filename().string();

                std::string suf = ".aln." + std::to_string(alnVersion);
                return (base_fn_2.empty()) ?
                        StringPair(base_fn_1 + suf, std::string()) :
                        StringPair(base_fn_1 + suf, base_fn_2 + suf);
        }

        /**
         * Get the fragment size
         * @return Fragment size
         */
        float getFragmentSize() const {
                return fragmentSize;
        }

        /**
         * Get the fragment standard deviation
         * @return Fragment standard deviation
         */
        float getFragmentStd() const {
                return fragmentStd;
        }

        /**
         * Filter uninformative alignments + estimated fragment size
         * @param dBG Const-reference to de Bruijn graph
         */
        void filterAlignments(const DBGraph& dBG);

        /**
         * Update read alignments taking to account graph manipulations
         * @param nodeTracker Collection of all graph manipulations
         * @param libary Read library
         */
        void parseAlignments(const NodeTracker& nodeTracker);

        /**
         * Count the number of reads that align to a node or edge
         * @param dBG Const-reference to de Bruijn graph
         * @param nodeCount Container to store node counts (output)
         * @param edgeCount Container to store edge counts (output)
         */
        void countAlignments(const DBGraph& dBG,
                             NodeMap<int>& nodeCount,
                             EdgeMap<int>& edgeCount);
};

// ============================================================================
// LIBRARY CONTAINER CLASS
// ============================================================================

class LibraryContainer : public std::vector<Library>
{
private:
        typedef std::pair<std::string, std::string> pairOfStrings;

public:
        /**
         * Constructor
         * @param manifestFilename Read filename
         */
        LibraryContainer(const std::string& manifestFilename);

        /**
         * Filter the informative alignments and estimate the insert size
         * @param dBG Const-reference to de Bruijn graph
         */
        void filterAlignments(const DBGraph& dBG) {
                for (Library& lib : *this)
                        lib.filterAlignments(dBG);
        }

        /**
         * Update read alignments taking to account graph manipulations
         * @param nodeTracker Collection of all graph manipulations
         */
        void parseAlignments(const NodeTracker& nodeTracker) {
                for (Library& lib : *this)
                        lib.parseAlignments(nodeTracker);
        }

        /**
         * Count the number of reads that align to a node or edge
         * @param dBG Const-reference to de Bruijn graph
         * @param nodeCount Container to store node counts (output)
         * @param edgeCount Container to store edge counts (output)
         */
        void countAlignments(const DBGraph& dBG,
                             NodeMap<int>& nodeCount,
                             EdgeMap<int>& edgeCount) {
                for (Library& lib : *this)
                        lib.countAlignments(dBG, nodeCount, edgeCount);
        }
};

#endif
