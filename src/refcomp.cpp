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
        FastAReader reader(fastaFilename);

        records.clear();
        FastARecord record;
        while (reader.getNextRecord(record))
                if (record.second.size() >= 250)
                        records.push_back(record);

        cout << "File: " << fastaFilename << "\n";
        /*for (const FastARecord& r : records) {
                const string& sequence = r.getSequence();
                cout << "\tSequence " << i << " length: " << sequence[i].size();

                size_t nonACTG = 0;
                for (size_t j = 0; j < sequence.size(); j++) {
                        char c = sequence[j];
                        if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
                                nonACTG++;
                }

                cout << " (" << nonACTG << " non-ACGT)" << endl;
        }*/
}

void RefComp::alignSequence(const string& seq,
                            vector<AlnSegment>& segment) const
{
        segment.clear();
        if (seq.size() < Kmer::getK())
                return;         // sequence too short

        // start by opening a CONTIG or BREAK
        NodePosPair prev = findNPP(seq.substr(0, Kmer::getK()));
        AlnSegmentType type = prev.isValid() ?
                AlnSegmentType::CONTIG : AlnSegmentType::BREAK;
        segment.push_back(AlnSegment(type));
        segment.back().setBegin(0, prev.getNodeID(), prev.getPosition());
        segment.back().addNode(prev.getNodeID());

        for (size_t i = 1; i < seq.size() + 1 - Kmer::getK(); i++) {
                NodePosPair curr = findNPP(seq.substr(i, Kmer::getK()));
                bool consecutive = dBG.consecutiveNPP(prev, curr);

                //if (curr.getNodeID() != prev.getNodeID())
                 //       cout << curr.getNodeID() << endl;

                if (type == AlnSegmentType::CONTIG && consecutive) {
                        // do we cross an arc
                        if ((prev.getPosition() + 1) != curr.getPosition())
                                segment.back().addNode(curr.getNodeID());
                }

                // NOTE: it is possible that both if-statements below are
                // executed within one loop iteration. Do not change the order!

                // encounter a breakpoint on a contig -> switch to break
                if (type == AlnSegmentType::CONTIG && !consecutive) {
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
                        segment.back().addNode(curr.getNodeID());
                }

                prev = curr;
        }

        // close open segment at the end of a sequence
        segment.back().setEnd(seq.size() + 1 - Kmer::getK(),
                              prev.getNodeID(), prev.getPosition()+1);
}

void RefComp::annotateBreakpoint(size_t i, const string& seq,
                                 vector<AlnSegment>& segment,
                                 GraphAligner& graphAln) const
{
        AlnSegment& ps = segment[i-1];
        AlnSegment& rs = segment[i];
        AlnSegment& ns = segment[i+1];

        if (rs.type != AlnSegmentType::BREAK)
                return;

        if (i == 0 || i == (segment.size() - 1)) {
                rs.type = AlnSegmentType::DELETION;
                return;
        }

        // go to the end of the previous CONTIG
        NodePosPair srcNPP = ps.getEnd();
        size_t srcSeq = ps.seqEnd;

        // rewind by at most MAX_NODES nodes or MAX_NUCLEOTIDES nucleotides
        size_t deltaBwd = 0, numRewind = 0;
        for (size_t j = 0; j < MAX_NODES_REWIND; j++) {
                // go to the end of the previous node
                size_t delta = srcNPP.getPosition();
                deltaBwd += delta;
                if (deltaBwd > MAX_NUCLEOTIDES_REWIND)
                        break;  // do not rewind more than MAX_NUCLEOTIDES
                if (deltaBwd >= ps.getSeqLength())
                        break;  // do not rewind entire sequence (or more)
                assert(j+1 < ps.nc.size());       // sanity check
                NodeID endNodeID = ps.nc[ps.nc.size() - 2 - j];
                SSNode endNode = dBG.getSSNode(endNodeID);
                srcNPP = NodePosPair(endNodeID, endNode.getMarginalLength());
                srcSeq -= delta;
                numRewind++;
        }

        // go to the start of the next CONTIG
        NodePosPair dstNPP = ns.getBegin();
        size_t dstSeq = ns.seqBegin;

        // forward by at most MAX_NODES nodes or MAX_NUCLEOTIDES nucleotides
        size_t deltaFwd = 0, numForward = 0;
        for (size_t j = 0; j < MAX_NODES_REWIND; j++) {
                // go to the start of the next node
                SSNode startNode = dBG.getSSNode(dstNPP.getNodeID());
                size_t delta = startNode.getMarginalLength() - dstNPP.getPosition();
                deltaFwd += delta;
                if (deltaFwd > MAX_NUCLEOTIDES_REWIND)
                        break;  // do not forward more than MAX_NUCLEOTIDES
                if (deltaFwd >= ns.getSeqLength())
                        break;  // do not forward entire sequence (or more)
                assert(j+1 < ns.nc.size());       // sanity check
                NodeID startNodeID = ns.nc[j+1];
                dstNPP = NodePosPair(startNodeID, 0);
                dstSeq += delta;
                numForward++;
        }

        // because the GraphAligner works with nucleotides rather than k-mers,
        // transform positions from k-mer space to nucleotide space
        NodePosPair srcNPP_n(srcNPP.getNodeID(), srcNPP.getPosition() + Kmer::getK() - 1);
        size_t srcSeq_n = srcSeq + Kmer::getK() - 1;
        NodePosPair dstNPP_n(dstNPP.getNodeID(), dstNPP.getPosition() + Kmer::getK() - 1);
        size_t dstSeq_n = dstSeq + Kmer::getK() - 1;

        // use the graph aligner to find the best path from srcNPP_n to dstNPP_n
        GraphAlignment ga;
        string gap = seq.substr(srcSeq_n, dstSeq_n - srcSeq_n);
        SearchRes res = graphAln.findPath(gap, srcNPP_n, dstNPP_n, ga);

        // no path found
        if (res != SearchRes::OPTIMAL || ga.getPath().empty())
                return;

        // adjust the end of the previous CONTIG
        for (size_t i = 0; i < numRewind; i++)
                ps.nc.pop_back();
        ps.setEnd(srcSeq, srcNPP.getNodeID(), srcNPP.getPosition());

        // adjust the PARALLEL path
        for (const auto& na : ga.getPath()) {
                if (na.nodeBegin == na.nodeEnd)
                        continue;       // skip empty aln (possible for first/last node)
                if (rs.nc.empty())      // set only the first time
                        rs.setBegin(srcSeq, na.nodeID, na.nodeBegin + 1 - Kmer::getK());
                rs.addNode(na.nodeID);
                rs.setEnd(dstSeq, na.nodeID, na.nodeEnd + 1 - Kmer::getK());
        }

        // adjust the beginning of the next CONTIG
        ns.setBegin(dstSeq, dstNPP.getNodeID(), dstNPP.getPosition());
        for (size_t i = 0; i < numForward; i++)
                ns.nc.pop_front();

        if (segment[i-1].seqEnd == segment[i+1].seqBegin) {
                rs.type = AlnSegmentType::INSERTION;
        } else {                   // ref is partially missing
                if (segment[i-1].getEnd() == segment[i+1].getBegin())
                        rs.type = AlnSegmentType::DELETION;
                else {
                        rs.type = AlnSegmentType::PARALLEL;
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
        const unsigned int& numThreads = 1; //settings.getNumThreads();
        //cout << "Number of threads: " << numThreads << endl;

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

void RefComp::alignSequences(vector< vector<AlnSegment> >& alignment)
{
        Util::startChrono();
        cout << "Aligning reference sequences to dBG... "; cout.flush();
        alignment.clear();
        alignment.resize(records.size());
        for (size_t i = 0; i < records.size(); i++)
                alignSequence(records[i].getSequence(), alignment[i]);
        cout << "done (" << Util::stopChronoStr() << ")" << endl;

        size_t totAlnSegments = 0;
        for (size_t i = 0; i < records.size(); i++)
                totAlnSegments += alignment[i].size();
        cout << "Aligned the reference " << records.size()
             << " sequence(s) to dBG using " << totAlnSegments << " segments\n";

        for (size_t i = 0; i < records.size(); i++) {
                stringstream output;
                annotateBreakpoints(records[i].getSequence(),
                                    alignment[i], output.str());

                output << "Non-contiguous alignment for sequence " << i << endl;
                printSegments(records[i].second, alignment[i]);
        }
}

void RefComp::printSegments(const string& seq,
                            const std::vector<AlnSegment>& alignment) const
{
        /*for (const AlnSegment& refSeg : alignment)
                cout << refSeg << endl;*/

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
        cout.precision(2);
        cout << "Sequence aligns with:\n" << fixed;
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

void RefComp::writeAlignedSeqs(const std::string& filename)
{
        // align sequences to the de Bruijn graph
        vector<vector<AlnSegment> > aln;
        alignSequences(aln);

        // write aligned contigs to disk
        ofstream ofs(filename.c_str());
        size_t contigID = 0;
        for (size_t i = 0; i < records.size(); i++) {
                for (const AlnSegment& refSeg : aln[i]) {
                        if (refSeg.type != AlnSegmentType::CONTIG)
                                continue;
                        ofs << ">contig_" << contigID++ << "\n";
                        size_t b = refSeg.seqBegin;
                        size_t e = refSeg.seqEnd;

                        Util::writeSeqWrap(ofs, records[i].second.substr(b, e-b), 60);
                }
        }
}

void RefComp::getReductions(vector<NodeChain>& red)
{
        // get the node chains
        vector<string> contID;
        vector<NodeChain> nodeChain;
        getTrueNodeChain(contID, nodeChain);

        /*ofstream ofs("nodechain.fasta");
        size_t contigID = 0;
        for (auto chain : nodeChain) {
                vector<NodeID> id;
                for (auto it : chain)
                        id.push_back(it);

                string str;
                dBG.convertNodesToString(id, str);

                ofs << ">contig_" << contigID++ << "\n";
                Util::writeSeqWrap(ofs, str, 60);
        }*/

        // get the node multiplicities
        NodeMap<int> nodeMult; EdgeMap<int> edgeMult;
        getTrueMultiplicity(nodeMult, edgeMult);

        size_t reductID = 0;
        for (const auto chain : nodeChain) {
                size_t startID = chain.size();
                for (size_t i = 0; i < chain.size(); i++) {
                        if (nodeMult[abs(chain[i])] != 1)
                                continue;

                        if (startID < i) {
                                red.push_back(NodeChain());
                                for (size_t j = startID; j <= i; j++)
                                        red.back().push_back(chain[j]);
                        }

                        startID = i;
                }
        }
}

void RefComp::getTrueNodeChain(vector<string>& contID,
                               vector<NodeChain>& nodeChain)
{
        // align the reference sequences to the DBG to obtain true node chains
        for (const FastARecord& r : records) {
                const string& id = r.getID();

                // FIXME: delete !!
                /*if (id != ">contig_9266" && id != ">contig_6866")
                        continue;*/

                const string& seq = r.second;
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
                                contID.push_back(id);
                                nodeChain.push_back(NodeChain());
                                nodeChain.back().setCount(1);
                                nodeChain.back().push_back(curr.getNodeID());
                        }

                        prev = curr;
                }
        }
}

NodeChain RefComp::getNodeChain(int seqID, size_t b, size_t e)
{
        const FastARecord& r = records[seqID];
        const string& id = r.getID();
        const string seq = r.second.substr(b, e-b);
        NodeChain result;

        // handle the other kmers
        NodePosPair prev;
        for (KmerIt it(seq); it.isValid(); it++) {
                Kmer kmer = it.getKmer();
                NodePosPair curr = table.find(kmer);

                if (!curr.isValid()) {
                        if (result.empty() || result.back() != 0)
                                result.push_back(0);
                        prev = curr;
                        continue;
                }

                if (dBG.consecutiveNPP(prev, curr)) {
                        if (curr.getNodeID() != prev.getNodeID())
                                result.push_back(curr.getNodeID());
                } else {
                        result.setCount(1);
                        result.push_back(curr.getNodeID());
                }

                prev = curr;
        }

        return result;
}

void RefComp::printSequence(int seqID, size_t b, size_t e)
{
        const FastARecord& r = records[seqID];
        const string& id = r.getID();
        const string seq = r.second.substr(b, e-b);
        cout << seq << endl;
}

NodeChain RefComp::getNodeChain(NodeID beginID, size_t len)
{
        vector<string> contID;
        vector<NodeChain> nodeChains;

        getTrueNodeChain(contID, nodeChains);

        NodeChain result;
        for (const auto& nc : nodeChains) {
                for (size_t i = 0; i < nc.size(); i++) {
                        if (nc[i] == beginID) {
                                size_t eo = min<size_t>(i + len, nc.size());
                                auto b = nc.begin() + i;
                                auto e = nc.begin() + eo;
                                result = NodeChain(b, e);
                        } if (nc[i] == -beginID) {
                                size_t bo = (i < len) ? 0 : i - len;
                                auto b = nc.begin() + bo;
                                auto e = nc.begin() + i + 1;
                                result = NodeChain(b, e);
                                result.revCompl();
                        }
                }
        }

        return result;
}

void RefComp::getTrueMultiplicity(NodeMap<int>& nodeMult,
                                  EdgeMap<int>& edgeMult)
{
        // create an empty node and edge map
        nodeMult = NodeMap<int>(dBG.getNumValidNodes());
        edgeMult = EdgeMap<int>(dBG.getNumValidArcs() / 2);

        for (size_t id = 1; id <= dBG.getNumNodes(); id++) {
                SSNode n = dBG.getSSNode(id);
                if (!n.isValid())
                        continue;
                nodeMult[id] = 0;

                for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++)
                        edgeMult[EdgeRep(it->getNodeID(), id)] = 0;
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++)
                        edgeMult[EdgeRep(id, it->getNodeID())] = 0;
        }

        // count the number of true kmers in each node
        for (const FastARecord& r : records) {
                const string& seq = r.second;
                // handle the first kmer separately
                KmerIt it(seq);
                NodePosPair prev = table.find(it.getKmer());
                if (prev.isValid())
                        nodeMult[prev.getNodeID()]++;

                // for all kmers in the reference sequence
                for (it++; it.isValid(); it++) {
                        Kmer kmer = it.getKmer();
                        NodePosPair curr = table.find(kmer);

                        // update kmer counters
                        if (curr.isValid())
                                nodeMult[curr.getNodeID()]++;

                        // get the arc multiplicity
                        if (dBG.crossesArc(prev, curr)) {
                                // arc from prev to curr
                                EdgeRep er(prev.getNodeID(), curr.getNodeID());
                                edgeMult[er]++;
                        }

                        prev = curr;
                }
        }

        size_t numTrueNodes = 0, numErrNodes = 0;

        // compute the multiplicity for each node and node statistics
        for (auto& it : nodeMult) {
                SSNode node = dBG.getSSNode(it.first);
                if (!node.isValid())
                        continue;

                double ML = node.getMarginalLength();
                it.second = round((double)it.second/ML);
                if (it.second > 0)
                        numTrueNodes++;
                else
                        numErrNodes++;
        }

        cout.precision(2);
        double perc = Util::toPercentage(numTrueNodes, nodeMult.size());
        cout << fixed << "Number of true nodes: " << numTrueNodes << "/"
             << nodeMult.size() << " (" << perc << "%)" << endl;
        perc = Util::toPercentage(numErrNodes, nodeMult.size());
        cout << fixed << "Number of false nodes: " << numErrNodes << "/"
             << nodeMult.size() << " (" << perc << "%)" << endl;

        // compute the edge statistics
        size_t numTrueArcs = 0, numErrArcs = 0;

        for (auto& it : edgeMult) {
                if (it.second > 0)
                        numTrueArcs++;
                else
                        numErrArcs++;
        }

        perc = Util::toPercentage(numTrueArcs, edgeMult.size());
        cout << fixed << "Number of true arcs: " << numTrueArcs << "/"
             << edgeMult.size() << " (" << perc << "%)" << endl;
        perc = Util::toPercentage(numErrArcs, edgeMult.size());
        cout << fixed << "Number of false arcs: " << numErrArcs << "/"
             << edgeMult.size() << " (" << perc << "%)" << endl;
}

void RefComp::getTrueMultiplicity(vector<vector<AlnSegment> >& aln,
                                  NodeMap<int>& nodeMult,
                                  EdgeMap<int>& edgeMult)
{
        // create an empty node and edge map
        nodeMult = NodeMap<int>(dBG.getNumValidNodes());
        edgeMult = EdgeMap<int>(dBG.getNumValidArcs() / 2);

        for (size_t id = 1; id <= dBG.getNumNodes(); id++) {
                SSNode n = dBG.getSSNode(id);
                if (!n.isValid())
                        continue;
                nodeMult[id] = 0;

                for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++)
                        edgeMult[EdgeRep(it->getNodeID(), id)] = 0;
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++)
                        edgeMult[EdgeRep(id, it->getNodeID())] = 0;
        }

        // count the number of true kmers in each node
        for (size_t k = 28; k < aln.size(); k++) {
        //for (const auto& a : aln) {
                const vector<AlnSegment> a = aln[k];
                for (size_t i = 0; i < a.size(); i++) {
                        const AlnSegment& s = a[i];

                        if ((s.type != AlnSegmentType::CONTIG) &&
                            (s.type != AlnSegmentType::PARALLEL))
                                continue;
                        bool continuation = ((i > 0) &&
                                (a[i-1].type == AlnSegmentType::CONTIG) ||
                                (a[i-1].type == AlnSegmentType::PARALLEL));

                        const deque<NodeID>& nc = s.nc;        // shorthand
                        for (size_t j = 0; j < nc.size(); j++) {
                                NodeID nodeID = nc[j];
                                SSNode n = dBG.getSSNode(nodeID);
                                int nl = n.getMarginalLength();
                                assert(n.isValid());

                                if (j == 0) {   // first node could be partial
                                        if (continuation) {
                                                EdgeRep er(a[i-1].nc.back(), nc[j]);
                                                if (dBG.edgeExists(er))
                                                        edgeMult[er]++;
                                        }
                                        nodeMult[nodeID] += (j+1 < nc.size()) ?
                                                nl - s.srcNodeBegin : s.dstNodeEnd - s.srcNodeBegin;
                                } else {
                                        EdgeRep er(nc[j-1], nc[j]);
                                        edgeMult[er]++;
                                        // last node could be partial
                                        nodeMult[nodeID] += (j+1 < nc.size()) ?
                                                nl : s.dstNodeEnd;
                                }
                        }
                }
        }

        size_t numTrueNodes = 0, numErrNodes = 0;

        // compute the multiplicity for each node and node statistics
        for (auto& it : nodeMult) {
                SSNode node = dBG.getSSNode(it.first);
                if (!node.isValid())
                        continue;

                double ML = node.getMarginalLength();
                it.second = round((double)it.second/ML);
                if (it.second > 0)
                        numTrueNodes++;
                else
                        numErrNodes++;
        }

        cout.precision(2);
        double perc = Util::toPercentage(numTrueNodes, nodeMult.size());
        cout << fixed << "Number of true nodes: " << numTrueNodes << "/"
             << nodeMult.size() << " (" << perc << "%)" << endl;
        perc = Util::toPercentage(numErrNodes, nodeMult.size());
        cout << fixed << "Number of false nodes: " << numErrNodes << "/"
             << nodeMult.size() << " (" << perc << "%)" << endl;

        // compute the edge statistics
        size_t numTrueArcs = 0, numErrArcs = 0;

        for (auto& it : edgeMult) {
                if (it.second > 0)
                        numTrueArcs++;
                else
                        numErrArcs++;
        }

        perc = Util::toPercentage(numTrueArcs, edgeMult.size());
        cout << fixed << "Number of true arcs: " << numTrueArcs << "/"
             << edgeMult.size() << " (" << perc << "%)" << endl;
        perc = Util::toPercentage(numErrArcs, edgeMult.size());
        cout << fixed << "Number of false arcs: " << numErrArcs << "/"
             << edgeMult.size() << " (" << perc << "%)" << endl;
}
