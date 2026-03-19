/***************************************************************************
 *   Copyright (C) 2014, 2015 Jan Fostier (jan.fostier@intec.ugent.be)     *
 *   This file is part of Brownie                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "global.h"
#include "kmeroverlaptable.h"
#include "kmeroverlap.h"
#include "settings.h"
#include "kmer/tstring.h"
#include "library.h"

#include "readfile/fasta.h"
#include "readfile/fastq.h"
#include "readfile/seqfile.h"

#include <deque>
#include <set>
#include <map>
#include <fstream>
#include <iostream>

using namespace std;

// ============================================================================
// KMER OVERLAP TABLE
// ============================================================================

KmerOverlapRef KmerOverlapTable::insert(const Kmer &kmer)
{
        // chose a representative kmer
        Kmer representative = kmer.getRepresentative();

        bool reverse = (kmer != representative);

        KmerOverlapPair val(representative, KmerOverlap());
        pair<KmerOverlapIt, bool> insResult = table.insert(val);
        KmerOverlapRef result(insResult.first, reverse);

        return result;
}

KmerOverlapRef KmerOverlapTable::find(const Kmer &kmer) const
{
        // chose a representative kmer
        Kmer representative = kmer.getRepresentative();

        bool reverse = (kmer != representative);
        return KmerOverlapRef(table.find(representative), reverse);
}

bool KmerOverlapTable::getLeftUniqueKmer(const KmerOverlapRef& rKmerRef,
                                  KmerOverlapRef& lKmerRef) const
{
        // initialise the right kmer reference to point to nothing
        lKmerRef = KmerOverlapRef(table.end(), false);

        char nucleotide;
        if (!rKmerRef.hasLeftUniqueOverlap(nucleotide))
                return false;

        Kmer lKmer = rKmerRef.getKmer();
        lKmer.pushNucleotideLeft(nucleotide);
        lKmerRef = find(lKmer);

        assert(lKmerRef.first != table.end());

        if (!lKmerRef.hasRightUniqueOverlap(nucleotide))
                return false;

        return true;
}

bool KmerOverlapTable::getRightUniqueKmer(const KmerOverlapRef& lKmerRef,
                                   KmerOverlapRef& rKmerRef) const
{
        // initialise the right kmer reference to point to nothing
        rKmerRef = KmerOverlapRef(table.end(), false);

        char nucleotide;
        if (!lKmerRef.hasRightUniqueOverlap(nucleotide))
                return false;

        Kmer rKmer = lKmerRef.getKmer();
        rKmer.pushNucleotideRight(nucleotide);
        rKmerRef = find(rKmer);

        assert(rKmerRef.first != table.end());

        if (!rKmerRef.hasLeftUniqueOverlap(nucleotide))
                return false;

        return true;
}

void KmerOverlapTable::convertKmersToString(const deque<KmerOverlapRef> &kmerSeq,
                                            string &output)
{
        output = kmerSeq[0].getKmer().str();

        for (size_t i = 1; i < kmerSeq.size(); i++)
                output.push_back(kmerSeq[i].getKmer().peekNucleotideRight());
}

void KmerOverlapTable::loadKmersFromDisc(const std::string& filename)
{
        // load all the kmers from disc
        ifstream ifs(filename.c_str(), ios::in | ios::binary);
        size_t numKmers;
        ifs.read((char*)&numKmers, sizeof(size_t));

        table.resize(numKmers);

        for (size_t i = 0; i < numKmers; i++) {
                Kmer kmer(ifs);
                insert(kmer);
        }

        ifs.close();
}

void KmerOverlapTable::parseRead(const FastQRecord& rr) const
{
        const string& read = rr.getRead();
        const string& qual = rr.getQual();

        // get out early
        if (read.size() < Kmer::getK())
                return;

        // reserve space for kmers and flags
        vector<KmerOverlapRef> refs(read.size() + 1 - Kmer::getK());

        // find the kmers in the table
        for (KmerIt it(read); it.isValid(); it++)
                refs[it.getOffset()] = find(it.getKmer());

        // now mark the overlap implied by the read
        //size_t lastIndex = 0;
        for (KmerIt it(read); it.isValid(); it++) {
                if (refs[it.getOffset()].first == table.end())
                        continue;
                if (it.hasRightOverlap())
                        if (refs[it.getOffset()+1].first != table.end())
                                refs[it.getOffset()].markRightOverlap(it.getRightOverlap());
                if (it.hasLeftOverlap())
                        if (refs[it.getOffset()-1].first != table.end())
                                refs[it.getOffset()].markLeftOverlap(it.getLeftOverlap());
        }
}

void KmerOverlapTable::workerThread(FastQReader& inputs)
{
        // local storage of reads
        vector<FastQRecord> readBuf;

        size_t chunkID;
        while (inputs.getNextChunk(readBuf, chunkID))
                for (const auto& readRecord : readBuf)
                        parseRead(readRecord);
}

void KmerOverlapTable::parseInputFiles(LibraryContainer &inputs)
{
        const unsigned int& numThreads = settings.getNumThreads();
        cout << "Number of threads: " << numThreads << endl;

        const size_t ws = settings.getThreadIOWorkSize();
        for (const Library& lib : inputs) {
                const auto [fn1, fn2] = lib.getFilename();

                FastQReader myReader(fn1, fn2);
                myReader.startReaderThread(ws, ws * settings.getNumThreads());

                // start worker threads
                vector<thread> workerThreads(numThreads);
                for (size_t i = 0; i < workerThreads.size(); i++)
                        workerThreads[i] = thread(&KmerOverlapTable::workerThread,
                                                  this, ref(myReader));

                // wait for worker threads to finish
                for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

                myReader.joinReaderThread();
        }
}

void KmerOverlapTable::extractNodes(const string& filename)
{
        ofstream ofsNode(filename + ".node.st1");
        ofstream ofsEdge(filename + ".edge.st1");

        string sequence;
        size_t numNodes = 0, numArcs = 0;
        for (KmerOverlapIt it = table.begin(); it != table.end(); it++) {

                KmerOverlapRef currKmer(it, false), nextKmer;

                if (numNodes % 1024 == 0)
                        (cout << "Extracting unitig " << numNodes
                              << " from graph\r").flush();

                // check if the node has been processed before
                if (currKmer.isProcessed()) continue;

                deque<KmerOverlapRef> kmerSeq;
                kmerSeq.push_back(currKmer);

                // extend node to the right
                while (getRightUniqueKmer(currKmer, nextKmer)) {

                        // check for a loop
                        if (nextKmer == kmerSeq.front())
                                break;
                        // check for a hairpin
                        if (nextKmer.first == kmerSeq.back().first)
                                break;

                        kmerSeq.push_back(nextKmer);
                        currKmer = nextKmer;
                }

                // extend node to the left
                currKmer = KmerOverlapRef(it, false);

                while (getLeftUniqueKmer(currKmer, nextKmer)) {

                        // check for a loop
                        if (nextKmer == kmerSeq.back())
                                break;
                        // check for a hairpin
                        if (nextKmer.first == kmerSeq.front().first)
                                break;

                        kmerSeq.push_front(nextKmer);
                        currKmer = nextKmer;
                }

                // mark all kmers as processed
                for (size_t i = 0; i < kmerSeq.size(); i++)
                        kmerSeq[i].setProcessed(true);

                convertKmersToString(kmerSeq, sequence);

                ofsNode << sequence << "\n";
                ofsEdge << (int)kmerSeq[0].getLeftOverlap() << "\t"
                        << (int)kmerSeq.back().getRightOverlap() << "\n";

                numNodes++;
                numArcs += kmerSeq[0].getNumLeftOverlap() +
                           kmerSeq.back().getNumRightOverlap();
        }

        ofsNode.close();
        ofsEdge.close();

        ofstream ofsMeta(filename + ".meta.st1");
        ofsMeta << Kmer::getK() << "\t" << numNodes << "\t" << numArcs << endl;
        ofsMeta.close();

        cout << "Extracted " << numNodes << " unitigs and "
             << numArcs << " arcs" << endl;
}
