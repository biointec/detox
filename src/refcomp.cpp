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

#include <iomanip>
#include <iostream>
#include <thread>

#include "refcomp.h"
#include "dbgraph.h"
#include "settings.h"
#include "graphaln.h"
#include "util.h"
#include "readfile/fastafile.h"

using namespace std;

std::ostream &operator<<(std::ostream &out, const AlnSegment &s)
{
        switch (s.type) {
                case AlnSegmentType::DELETION:
                        out << "DELETION  ";
                        break;
                case AlnSegmentType::INSERTION:
                        out << "INSERTION ";
                        break;
                case AlnSegmentType::PARALLEL:
                        out << "PARALLEL  ";
                        break;
                case AlnSegmentType::BREAK:
                        out << "BREAK     ";
                        break;
                case AlnSegmentType::CONTIG:
                        out << "CONTIG    ";
                        break;
        }

        out << "[" <<  s.seqBegin << " - " << s.seqEnd << "[ ";
        if (s.srcNodeID == 0)
                out << "**" << " - ";
        else
                out << s.srcNodeID << " (" << s.srcNodeBegin << ") - ";

        if (s.dstNodeID == 0)
                out << "**";
        else
                out << s.dstNodeID << " (" << s.dstNodeEnd << ")";

        return out;
}

// ============================================================================
// REFERENCE COMPARISON CLASS
// ============================================================================

void RefComp::readFastaFile(const string& fastaFilename)
{
        sequence.clear();

        FastAFile ifs(false);
        ifs.open(fastaFilename.c_str());

        string read;
        while (ifs.getNextRead(read))
                sequence.push_back(read);

        cout << "File: " << fastaFilename << "\n";
        for (size_t i = 0; i < sequence.size(); i++) {
                cout << "\tSequence " << i << " length: "
                     << sequence[i].size();

                size_t nonACTG = 0;
                for (size_t j = 0; j < sequence[i].size(); j++) {
                        char c = sequence[i][j];
                        if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
                                nonACTG++;
                }

                cout << " (" << nonACTG << " non-ACGT)" << endl;
        }

        ifs.close();
}

void RefComp::alignSequence(const string& seq,
                            vector<AlnSegment>& segment) const
{
        /*segment.clear();
        if (seq.size() < Kmer::getK())
                return;         // sequence too short

        // start by opening a CONTIG or BREAK
        NodePosPair prev = findNPP(seq.substr(0, Kmer::getK()));
        AlnSegmentType type = prev.isValid() ?
                AlnSegmentType::CONTIG : AlnSegmentType::BREAK;
        segment.push_back(AlnSegment(type));
        segment.back().setBegin(0, prev.getNodeID(), prev.getPosition());

        for (size_t i = 1; i < seq.size() + 1 - Kmer::getK(); i++) {
                NodePosPair curr = findNPP(seq.substr(i, Kmer::getK()));

                // NOTE: it is possible that both if-statements below are
                // executed within one loop iteration. Do not change the order!

                // encounter a breakpoint on a contig -> switch to break
                if (type == AlnSegmentType::CONTIG && !dBG.consecutiveNPP(prev, curr)) {
                        // close the contig
                        segment.back().setEnd(i, prev.getNodeID(),
                                              prev.getPosition()+1);

                        // push a breakpoint
                        type = AlnSegmentType::BREAK;
                        segment.push_back(AlnSegment(type));
                        segment.back().setBegin(i, 0, 0);
                }

                // encounter a valid kmer with no open segment -> switch to contig
                if (type == AlnSegmentType::BREAK && curr.isValid()) {
                        // close the breakpoint
                        segment.back().setEnd(i, 0, 0);

                        // push a contig
                        type = AlnSegmentType::CONTIG;
                        segment.push_back(AlnSegment(type));
                        segment.back().setBegin(i, curr.getNodeID(),
                                                curr.getPosition());
                }

                prev = curr;
        }

        // close open segment at the end of a sequence
        segment.back().setEnd(seq.size() + 1 - Kmer::getK(),
                              prev.getNodeID(), prev.getPosition()+1);*/
}

void RefComp::annotateBreakpoint(size_t i, const string& seq,
                                 vector<AlnSegment>& segment,
                                 GraphAligner& graphAln) const
{
        AlnSegment& rs = segment[i];

        if (rs.type != AlnSegmentType::BREAK)
                return;

        if (i == 0 || i == (segment.size() - 1)) {
                rs.type = AlnSegmentType::DELETION;
                return;
        }

        NodePosPair srcNPP = segment[i-1].getEnd();
        size_t srcSeq = segment[i-1].seqEnd;

        NodePosPair dstNPP = segment[i+1].getBegin();
        size_t dstSeq = segment[i+1].seqBegin;

        if (dstNPP.getPosition() < Kmer::getK() - 1) {
                size_t delta = Kmer::getK() - 1 - dstNPP.getPosition();
                dstNPP.setPosition(Kmer::getK() - 1);
                dstSeq += delta;
        }

        GraphAlignment ga;
        string gap = seq.substr(srcSeq, dstSeq - srcSeq);
        SearchRes res = graphAln.findPath(gap, srcNPP, dstNPP, ga);

        if (res != SearchRes::OPTIMAL || ga.getPath().empty())
                return;

        if (segment[i-1].seqEnd == segment[i+1].seqBegin) {
                rs.type = AlnSegmentType::INSERTION;
        } else {                   // ref is partially missing
                if (segment[i-1].getEnd() == segment[i+1].getBegin())
                        rs.type = AlnSegmentType::DELETION;
                else {
                        rs.type = AlnSegmentType::PARALLEL;
                        //cout << seq.substr(srcSeq, dstSeq - srcSeq) << endl;
                }
        }
}

void RefComp::annotateBreakpointThread(size_t myID, WorkLoadBalancer& wlb,
                                       const string& seq,
                                       vector<AlnSegment>& segment) const
{
        // thread-local auxiliary variables
        GraphAligner graphAln(dBG);

        size_t chunkBegin, chunkEnd;
        while (wlb.getChunk(chunkBegin, chunkEnd)) {
                for (size_t i = chunkBegin; i < chunkEnd; i++)
                        annotateBreakpoint(i, seq, segment, graphAln);
        }
}

void RefComp::annotateBreakpoints(const string& seq,
                                  vector<AlnSegment>& segment,
                                  const string& output) const
{
        const unsigned int& numThreads = settings.getNumThreads();
        cout << "Number of threads: " << numThreads << endl;

        size_t targetWL = max<size_t>(1ull, segment.size() / 10 / numThreads);
        WorkLoadBalancer wlb(0, segment.size(), targetWL, output);

        // start worker threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&RefComp::annotateBreakpointThread,
                                          this, i, ref(wlb), cref(seq), ref(segment));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));
}

void RefComp::alignSequences()
{
        Util::startChrono();
        cout << "Aligning reference sequences to dBG... "; cout.flush();
        vector< vector<AlnSegment> > alignment(sequence.size());
        for (size_t i = 0; i < sequence.size(); i++)
                alignSequence(sequence[i], alignment[i]);
        cout << "done (" << Util::stopChronoStr() << ")" << endl;

        size_t totAlnSegments = 0;
        for (size_t i = 0; i < sequence.size(); i++)
                totAlnSegments += alignment[i].size();
        cout << "Aligned the reference sequence(s) to dBG using "
             << totAlnSegments << " segments." << endl;

        for (size_t i = 0; i < sequence.size(); i++) {
                Util::startChrono();
                stringstream output;
                output << "Annotating alignment for sequence " << i << "... ";
                annotateBreakpoints(sequence[i], alignment[i], output.str());
                cout << "done (" << Util::stopChronoStr() << ")" << endl;
                printSegments(sequence[i], alignment[i]);
        }
}

void RefComp::printSegments(const string& seq,
                            const std::vector<AlnSegment>& alignment) const
{
        // for (const AlnSegment& refSeg : alignment)
        //        cout << refSeg << endl;

        size_t numContigs = 0, numParallel = 0, numInserts = 0, numDeletions = 0, numBreaks = 0;
        size_t lenContigs = 0, lenParallel = 0, lenInserts = 0, lenDeletions = 0, lenBreaks = 0;
        for (const AlnSegment& refSeg : alignment) {
                switch (refSeg.type) {
                        case AlnSegmentType::CONTIG:
                                numContigs++;
                                lenContigs += refSeg.getSeqLength();
                                break;
                        case AlnSegmentType::PARALLEL:
                                numParallel++;
                                lenParallel += refSeg.getSeqLength();
                                break;
                        case AlnSegmentType::INSERTION:
                                numInserts++;
                                lenInserts += refSeg.getSeqLength();
                                break;
                        case AlnSegmentType::DELETION:
                                numDeletions++;
                                lenDeletions += refSeg.getSeqLength();
                                break;
                        case AlnSegmentType::BREAK:
                                numBreaks++;
                                lenBreaks += refSeg.getSeqLength();
                                break;
                }
        }

        size_t margLen = seq.length() + 1 - Kmer::getK();
        cout << "Sequence aligns with:\n";
        cout << "\tContigs: "   << numContigs   << ", total length: " << lenContigs
             << " (" << Util::toPercentage(lenContigs, margLen) << "%)" << endl;
        cout << "\tParallel: "  << numParallel  << ", total length: " << lenParallel
             << " (" << Util::toPercentage(lenParallel, margLen) << "%)" << endl;
        cout << "\tInserts: "   << numInserts   << ", total length: " << lenInserts
             << " (" << Util::toPercentage(lenInserts, margLen) << "%)" << endl;
        cout << "\tDeletions: " << numDeletions << ", total length: " << lenDeletions
             << " (" << Util::toPercentage(lenDeletions, margLen) << "%)" << endl;
        cout << "\tBreaks: "    << numBreaks    << ", total length: " << lenBreaks
             << " (" << Util::toPercentage(lenBreaks, margLen) << "%)" << endl;
}

void RefComp::getTrueNodeChain(vector<NodeChain>& nodeChain)
{
        // align the reference sequences to the DBG to figure out true node chains
        /*for (const string& seq : sequence) {
                // handle the other kmers
                NodePosPair prev;
                for (KmerIt it(seq); it.isValid(); it++) {
                        Kmer kmer = it.getKmer();
                        NodePosPair curr = table.find(kmer);

                        if (!curr.isValid()) {
                                prev = curr;
                                continue;
                        }

                        if (dBG.consecutiveNPP(prev, curr)) {
                                if (curr.getNodeID() != prev.getNodeID())
                                        nodeChain.back().push_back(curr.getNodeID());
                        } else {
                                nodeChain.push_back(NodeChain());
                                nodeChain.back().setCount(1);
                                nodeChain.back().push_back(curr.getNodeID());
                        }

                        prev = curr;
                }
        }*/
}

void RefComp::getTrueMultiplicity(const KmerNPPTable& table,
                                  vector<int>& nodeMult,
                                  vector<int>& edgeMult)
{
        nodeMult.clear();
        nodeMult.resize(dBG.getNumNodes()+1, 0);
        edgeMult.clear();
        edgeMult.resize(dBG.getNumArcs()+1, 0);

        // count the number of true kmers in each node
        for (const string& seq : sequence) {
                // handle the first kmer separately
                KmerIt it(seq);
                NodePosPair prev = table.find(it.getKmer());
                if (prev.isValid())
                        nodeMult[abs(prev.getNodeID())]++;

                // for all kmers in the reference sequence
                for (it++; it.isValid(); it++) {
                        Kmer kmer = it.getKmer();
                        NodePosPair curr = table.find(kmer);

                        // update kmer counters
                        if (curr.isValid())
                                nodeMult[abs(curr.getNodeID())]++;

                        // get the arc multiplicity
                        if (dBG.crossesArc(prev, curr)) {
                                // arc from prev to curr
                                SSNode left = dBG.getSSNode(prev.getNodeID());
                                Arc* arc = left.rightArc(curr.getNodeID());
                                edgeMult[dBG.getArcID(arc)]++;

                                // arc from curr to prev (don't count palindromic arcs twice!)
                                if (prev.getNodeID() != -curr.getNodeID()) {
                                        SSNode right = dBG.getSSNode(curr.getNodeID());
                                        arc = right.leftArc(prev.getNodeID());
                                        edgeMult[dBG.getArcID(arc)]++;
                                }
                        }

                        prev = curr;
                }
        }

        size_t numTrueNodes = 0, numErrNodes = 0;

        // compute the multiplicity for each node
        for (size_t i = 1; i <= dBG.getNumNodes(); i++) {
                SSNode node = dBG.getSSNode(i);
                if (!node.isValid())
                        continue;

                double ML = node.getMarginalLength();
                nodeMult[i] = round((double)nodeMult[i]/ML);
                if (nodeMult[i] > 0)
                        numTrueNodes++;
                else
                        numErrNodes++;
        }

        double perc = Util::toPercentage(numTrueNodes, dBG.getNumValidNodes());
        cout << "Number of true nodes: " << numTrueNodes << "/"
             << dBG.getNumValidNodes() << " (" << perc << "%)" << endl;
        perc = Util::toPercentage(numErrNodes, dBG.getNumValidNodes());
        cout << "Number of false nodes: " << numErrNodes << "/"
             << dBG.getNumValidNodes() << " (" << perc << "%)" << endl;

        size_t numTrueArcs = 0, numErrArcs = 0;

        // compute the multiplicity for each node
        for (size_t i = 1; i <= dBG.getNumArcs(); i++) {
                Arc& arc = dBG.getArc(i);
                if (!arc.isValid())
                        continue;

                if (edgeMult[i] > 0)
                        numTrueArcs++;
                else
                        numErrArcs++;
        }

        perc = Util::toPercentage(numTrueArcs, dBG.getNumValidArcs());
        cout << "Number of true arcs: " << numTrueArcs << "/"
             << dBG.getNumValidArcs() << " (" << perc << "%)" << endl;
        perc = Util::toPercentage(numErrArcs, dBG.getNumValidArcs());
        cout << "Number of false arcs: " << numErrArcs << "/"
             << dBG.getNumValidArcs() << " (" << perc << "%)" << endl;
}
