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

#include "readfile/fastq.h"
#include "nodechain.h"
#include "library.h"
#include "dbgraph.h"
#include "nodetracker.h"

using namespace std;

// ============================================================================
// LIBRARY CLASS
// ============================================================================

ReadOrientation Library::inferReadOrientation(const DBGraph& dBG,
                                              size_t numReads) const
{
        size_t ff = 0, rf = 0, fr = 0;

        ifstream ifs1(filename.first);
        ifstream ifs2(filename.second);

        string line1, line2; size_t counter = 0;
        while (getline(ifs1, line1) && getline(ifs2, line2)) {
                int begR1, endR1, begAln1, endAln1, RL1;
                int begR2, endR2, begAln2, endAln2, RL2;
                NodeChain R1, R2;

                istringstream iss(line1);
                NodeID nodeID;
                iss >> begR1 >> endR1 >> begAln1 >> endAln1 >> RL1;
                while (iss >> nodeID)
                        R1.push_back(nodeID);

                iss = istringstream(line2);
                iss >> begR2 >> endR2 >> begAln2 >> endAln2 >> RL2;
                while (iss >> nodeID)
                        R2.push_back(nodeID);

                if (R1.size() != 1 || R2.size() != 1)
                        continue;

                const NodeID& id1 = R1[0]; const NodeID& id2 = R2[0];
                if (NodeRep(id1) != NodeRep(id2))
                        continue;

                if (id1 == id2) {       // ff case
                        ff++;
                } else {                // rf or fr case
                        // change the orientation of R2
                        // if it comes after R1, orientation is fr, else rf
                        int nl = dBG.getSSNode(id1).length();
                        if (begAln1 < (nl - begAln2))
                                fr++;
                        else
                                rf++;
                }

                if (++counter == numReads)
                        break;
        }

        if (ff > rf)    // ff or fr
                return (ff > fr) ? FF : FR;
        else            // rf or fr
                return (rf > fr) ? RF : FR;
}

void Library::filterAlignments(const DBGraph& dBG)
{
        ReadOrientation ro = inferReadOrientation(dBG);
        if (ro == FF)
                cout << "Paired-end read orientation: FF --->   --->\n";
        if (ro == RF)
                cout << "Paired-end read orientation: RF <---   --->\n";
        if (ro == FR)
                cout << "Paired-end read orientation: FR --->   <---\n";

        const auto [fn1In, fn2In] = getAlnFilename();
        alnVersion++;
        const auto [fn1Out, fn2Out] = getAlnFilename();

        ifstream ifs1(fn1In);
        ifstream ifs2(fn2In);

        ofstream ofs1(fn1Out);
        ofstream ofs2(fn2Out);

        MeanStd fragSizeEst;
        size_t readsIn = 0, readsOut = 0;

        string line1, line2; int minSpace;
        while (getline(ifs1, line1) && getline(ifs2, line2)) {
                // read in both reads
                int begR1, endR1, begAln1, endAln1, RL1;
                int begR2, endR2, begAln2, endAln2, RL2;
                NodeChain R1, R2;

                istringstream iss(line1);
                NodeID nodeID;
                iss >> begR1 >> endR1 >> begAln1 >> endAln1 >> RL1;
                while (iss >> nodeID)
                        R1.push_back(nodeID);

                iss = istringstream(line2);
                iss >> begR2 >> endR2 >> begAln2 >> endAln2 >> RL2;
                while (iss >> nodeID)
                        R2.push_back(nodeID);

                // update sigma once in a while
                if ((readsIn++ % 4096) == 0) {
                        minSpace = (fragSizeEst.getNumObs() < 10) ?
                        0 : int(fragSizeEst.getMean() + 3.0 * fragSizeEst.getStd());
                }

                // orient reads as FF
                if (ro == RF)
                        dBG.reverseRead(R1, begR1, endR1, begAln1, endAln1, RL1);
                else if (ro == FR)
                        dBG.reverseRead(R2, begR2, endR2, begAln2, endAln2, RL2);

                // skip uninformative reads
                if (R1.empty() && R2.size() <= 1)
                        continue;
                if (R2.empty() && R1.size() <= 1)
                        continue;

                // estimate the fragment length
                if ((R1.size() == 1) && (R1 == R2)) {
                        int b = begAln1 - begR1;
                        int e = endAln2 + RL2 - endR2;
                        int fragSize = e - b;

                        int minNL = b + minSpace;

                        if ((dBG.getSSNode(R1[0]).length() >= minNL) && (fragSize > 0))
                                fragSizeEst.addObservation(fragSize);

                        continue;       // skip uninformative read
                }

                ofs1 << begAln1 << "\t" << endAln1;
                for (size_t i = 0; i < R1.size(); i++)
                        ofs1 << "\t" << R1[i];
                ofs1 << "\n";
                ofs2 << begAln2 << "\t" << endAln2;
                for (size_t i = 0; i < R2.size(); i++)
                        ofs2 << "\t" << R2[i];
                ofs2 << "\n";

                readsOut++;
        }

        fragmentSize = fragSizeEst.getMean();
        fragmentStd = fragSizeEst.getStd();

        cout << "Estimated fragment size: " << std::fixed
             << fragmentSize
             << " (std = " << fragmentStd << "), based on "
             << fragSizeEst.getNumObs() << " observations" << endl;

        cout << "Parsed " << readsIn << " read alignments. Kept "
             << readsOut << " informative alignments (version 0)\n";
}

void Library::parseAlignments(const NodeTracker& nodeTracker)
{
        size_t readsIn = 0, readsOut = 0;
        const auto [fn1In, fn2In] = getAlnFilename();
        alnVersion++;
        const auto [fn1Out, fn2Out] = getAlnFilename();

        ifstream ifs1(fn1In);
        ifstream ifs2(fn2In);

        ofstream ofs1(fn1Out);
        ofstream ofs2(fn2Out);

        string line1, line2;
        while (getline(ifs1, line1) && getline(ifs2, line2)) {
                readsIn++;
                int begAln1, endAln1, begAln2, endAln2;
                NodeChain R1, R2;

                istringstream iss(line1);
                NodeID nodeID;
                iss >> begAln1 >> endAln1;
                while (iss >> nodeID)
                        R1.push_back(nodeID);

                iss = istringstream(line2);
                iss >> begAln2 >> endAln2;
                while (iss >> nodeID)
                        R2.push_back(nodeID);

                nodeTracker.updateNodeChain(R1, begAln1, endAln1);
                nodeTracker.updateNodeChain(R2, begAln2, endAln2);

                if (R1.empty() && R2.size() <= 1)
                        continue;
                if (R2.empty() && R1.size() <= 1)
                        continue;
                if (R1.size() == 1 && R1 == R2)
                        continue;

                ofs1 << begAln1 << "\t" << endAln1;
                for (size_t i = 0; i < R1.size(); i++)
                        ofs1 << "\t" << R1[i];
                ofs1 << "\n";
                ofs2 << begAln2 << "\t" << endAln2;
                for (size_t i = 0; i < R2.size(); i++)
                        ofs2 << "\t" << R2[i];
                ofs2 << "\n";

                readsOut++;
        }

        cout << "Parsed " << readsIn << " read alignments. Kept "
             << readsOut << " informative alignments (version "
             << alnVersion << ")\n";
}

void Library::countAlignments(const DBGraph& dBG,
                              NodeMap<int>& nodeCount,
                              EdgeMap<int>& edgeCount)
{
        const auto [fn1In, fn2In] = getAlnFilename();

        ifstream ifs1(fn1In);
        ifstream ifs2(fn2In);

        string line1, line2;
        while (getline(ifs1, line1) && getline(ifs2, line2)) {
                int begAln1, endAln1, begAln2, endAln2;
                NodeChain R1, R2;

                istringstream iss(line1);
                NodeID nodeID;
                iss >> begAln1 >> endAln1;
                while (iss >> nodeID)
                        R1.push_back(nodeID);

                iss = istringstream(line2);
                iss >> begAln2 >> endAln2;
                while (iss >> nodeID)
                        R2.push_back(nodeID);

                for (NodeID id : R1)
                        nodeCount[id]++;
                for (size_t i = 1; i < R1.size(); i++)
                        edgeCount[EdgeRep(R1[i-1], R1[i])]++;

                for (NodeID id : R2)
                        nodeCount[id]++;
                for (size_t i = 1; i < R2.size(); i++)
                        edgeCount[EdgeRep(R2[i-1], R2[i])]++;
        }
}

// ============================================================================
// LIBRARY CONTAINER CLASS
// ============================================================================

LibraryContainer::LibraryContainer(const std::string& mf)
{
        try {                           // option A: mf is a FastQ file
                FastQReader reader(mf, "");
                push_back(Library(reader.getFilename(),
                                  reader.getBaseFilename()));
        } catch (exception& e) {        // option B: mf is a manifest file
                ifstream ifs(mf);
                if (!ifs)
                        throw runtime_error("cannot open file " + mf);

                string line;
                while (getline(ifs, line)) {
                        if (line.empty())
                                continue;
                        const string tokens = " \t,;";
                        size_t first = line.find_first_of(tokens);
                        size_t last  = line.find_last_of(tokens);

                        // no delimiter -> single-ended reads
                        if (last >= line.size()) {
                                FastQReader reader(line, "");
                                push_back(Library(reader.getFilename(),
                                                  reader.getBaseFilename()));
                        } else {        // paired-end reads
                                string file1 = line.substr(0, first);
                                string file2 = line.substr(last + 1);
                                FastQReader reader(file1, file2);
                                push_back(Library(reader.getFilename(),
                                                  reader.getBaseFilename()));
                        }
                }
        }
}
