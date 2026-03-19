/******************************************************************************
 *   Copyright (C) 2014 - 2026 Jan Fostier (jan.fostier@ugent.be)             *
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

#include <cstdlib>
#include <numeric>

#include "dbgraph.h"
#include "settings.h"
#include "refcomp.h"
#include "crfmult.h"
#include "readaln.h"
#include "nodetracker.h"
#include "kmertable.h"
#include "kmeroverlap.h"
#include "kmeroverlaptable.h"

#include "pathinfo.h"

#include "nninference.h"

using namespace std;

void populateNodeMult(NodeMap<Multiplicity>& nodeMult,
                      const vector<NodeRep>& nodes)
{
        nodeMult = NodeMap<Multiplicity>(nodes.size());
        for (size_t i = 0; i < nodes.size(); i++)
                nodeMult[nodes[i]] = Multiplicity();
}

void populateEdgeMult(EdgeMap<Multiplicity>& edgeMult,
                      const vector<EdgeRep>& edges)
{
        edgeMult = EdgeMap<Multiplicity>(edges.size());
        for (size_t i = 0; i < edges.size(); i++)
                edgeMult[edges[i]] = Multiplicity();
}

template<class T>       // in practice T = NodeMap<int> or EdgeMap<int>
void computeAccuracy(const T& trueMult, const T& estMult)
{
        size_t wrong = 0, total = 0;
        size_t zeroZero = 0, zeroOne = 0, zeroMore = 0,
               oneZero = 0, oneOne = 0, oneMore = 0;
        for (auto [id, tm] : trueMult)
        {
                int em = get(estMult, id);

                if (tm != em && (tm < 5 || em < 5))
                        wrong++;
                if (tm == 0 && em == 0)
                        zeroZero++;
                if (tm == 0 && em == 1)
                        zeroOne++;
                if (tm == 0 && em > 1)
                        zeroMore++;
                if (tm == 1 && em == 0)
                        oneZero++;
                if (tm == 1 && em == 1)
                        oneOne++;
                if (tm == 1 && em > 1)
                        oneMore++;

                total++;
        }
        cout << "Incorrect: " << wrong << "/" << total << " ("
             << 100.0*double(wrong)/double(total) << "\%)" << endl;
        cout << "   est| 0\t1\t2+\n";
        cout << "---------------------------\n";
        cout << "true 0| " << zeroZero << "\t" << zeroOne << "\t" << zeroMore << "\n"
             << "true 1| " << oneZero << "\t" << oneOne << "\t" << oneMore << endl;
}

void testAccuracy(DBGraph& dBG, Settings& settings, LibraryContainer& library)
{
        // load node/edge coverage model
        CovModel nodeModel(settings.getStage3NodeModelFilename());
        CovModel edgeModel(settings.getStage3EdgeModelFilename());

        cout << "Node model: " << nodeModel << endl;
        cout << "Edge model: " << edgeModel << endl;

        KmerNPPTable table(dBG, true);
        dBG.createKmerNPPTable(table);
        RefComp refComp(dBG, settings, table, "genome.fasta");
        NodeMap<int> trueNodeMult; EdgeMap<int> trueEdgeMult;
        refComp.getTrueMultiplicity(trueNodeMult, trueEdgeMult);

        NodeMap<Multiplicity> nodeMult; EdgeMap<Multiplicity> edgeMult;
        populateNodeMult(nodeMult, dBG.getNodeReps(dBG.getNumValidNodes()));
        populateEdgeMult(edgeMult, dBG.getEdgeReps(dBG.getNumValidArcs()));

        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMargin(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.getNumThreads(),
                          settings.getCRFWorkSize());

        myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel);

        ofstream ofs("truemult.node");
        for (auto [id, tm] : trueNodeMult) {
                SSNode n = dBG.getSSNode(id);
                string s = n.getSequence();
                int nCG = 0;
                for (int i = 0; i < s.size(); i++) {
                        if (s[i] == 'C' || s[i] == 'G')
                                nCG++;
                }
                float CGperc = (float)nCG / (float)s.size();

                ofs << id << "\t"
                    << n.getAvgCov() << "\t"
                    << n.getMarginalLength() << "\t"
                    << CGperc << "\t"
                    << tm << "\n";
        }
        ofs.close();

        ofs.open("truemult.edge");
        for (auto [id, em] : trueEdgeMult) {
                Arc& a = dBG.getArc(id);
                SSNode src = dBG.getSSNode(id.getSrcID());
                string s = src.getRightKmer().str();
                SSNode dst = dBG.getSSNode(id.getDstID());
                char c = dst.peekNucleotideMarginalLeft();
                s += c;
                int nCG = 0;
                for (int i = 0; i < s.size(); i++) {
                        if (s[i] == 'C' || s[i] == 'G')
                                nCG++;
                }
                float CGperc = (float)nCG / (float)s.size();

                ofs << id.getSrcID() << "\t" << id.getDstID() << "\t"
                    << a.getCov() << "\t"
                    << CGperc << "\t"
                    << em << "\n";
        }
        ofs.close();

        size_t wrong = 0, total = 0;
        size_t zeroZero = 0, zeroOne = 0, zeroMore = 0,
               oneZero = 0, oneOne = 0, oneMore = 0;
        for (auto [id, tm] : trueNodeMult) {
                /*if (tm > 1)
                        continue;*/

                int em = nodeMult[id].getExpMult();

                SSNode n = dBG.getSSNode(id);
                /*if ((tm != em) && (n.getMarginalLength() > 250))
                        cout << "Node " << id << " (" << n.getMarginalLength() << "), est: " << em << ", true: " << tm << endl;*/

                if (tm != em && (tm < 5 || em < 5))
                        wrong++;
                if (tm == 0 && em == 0)
                        zeroZero++;
                if (tm == 0 && em == 1)
                        zeroOne++;
                if (tm == 0 && em > 1)
                        zeroMore++;
                if (tm == 1 && em == 0)
                        oneZero++;
                if (tm == 1 && em == 1)
                        oneOne++;
                if (tm == 1 && em > 1)
                        oneMore++;

                total++;
        }
        cout << "Incorrect: " << wrong << "/" << total << " ("
             << 100.0*double(wrong)/double(total) << "\%)" << endl;
        cout << zeroZero << "\t" << zeroOne << "\t" << zeroMore << "\n"
             << oneZero << "\t" << oneOne << "\t" << oneMore << endl;
}

void trainModel(DBGraph& dBG, Settings& settings,
                CovModel& nodeModel, CovModel& edgeModel, bool trainNodeODF)
{
        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMargin(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.getNumThreads(),
                          settings.getCRFWorkSize());

        NodeMap<Multiplicity> nodeMult;
        EdgeMap<Multiplicity> edgeMult;
        size_t trainSize = settings.getEMTrainingSize();
        populateNodeMult(nodeMult, dBG.getNodeReps(trainSize));
        populateEdgeMult(edgeMult, dBG.getEdgeReps(trainSize));

        int nIters = myCRFMult.computeMultEM(nodeMult, nodeModel,
                                             edgeMult, edgeModel,
                                             trainNodeODF,
                                             settings.getEMConvEps(),
                                             settings.getEMMaxIter());

        if (nIters <= settings.getEMMaxIter())
                cout << "EM algorithm converged after "
                     << nIters << " iterations\n";
        else
                cout << "WARNING: maximum number of iterations reached. "
                     << "Convergence is not guaranteed." << endl;

        map<int, double> nodeHist;
        for (const auto& it : nodeMult) {
                SSNode node = dBG.getSSNode(it.first);
                double f = node.getAvgCov() - floor(node.getAvgCov());
                nodeHist[node.getAvgCov()] += (1.0-f);
                nodeHist[node.getAvgCov() + 1] += f;
        }

        nodeModel.writeGnuplot("histogram.node", nodeHist);

        map<int, double> edgeHist;
        for (const auto& it : edgeMult) {
                Arc& arc = dBG.getArc(dBG.getArcID(it.first));
                double f = arc.getCov() - floor(arc.getCov());
                edgeHist[arc.getCov()] += (1.0-f);
                edgeHist[arc.getCov() + 1] += f;
        }

        edgeModel.writeGnuplot("histogram.edge", edgeHist);
}

void computeMultiplicities(DBGraph& dBG, Settings& settings,
                           NodeMap<Multiplicity>& nodeMult,
                           EdgeMap<Multiplicity>& edgeMult)
{
        // get the node multiplicities
        CovModel nodeModel(settings.getStage3NodeModelFilename());
        CovModel edgeModel(settings.getStage3EdgeModelFilename());

        cout << "Node model: " << nodeModel << endl;
        cout << "Edge model: " << edgeModel << endl;

        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMargin(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.getNumThreads(),
                          settings.getCRFWorkSize());

        myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel);
}

void computeAllMultiplicities(DBGraph& dBG, Settings& settings,
                              NodeMap<Multiplicity>& nodeMult,
                              EdgeMap<Multiplicity>& edgeMult)
{
        populateNodeMult(nodeMult, dBG.getNodeReps(dBG.getNumValidNodes()));
        populateEdgeMult(edgeMult, dBG.getEdgeReps(dBG.getNumValidArcs()));
        computeMultiplicities(dBG, settings, nodeMult, edgeMult);
}

void stageOne(Settings& settings, LibraryContainer& libraries)
{
        cout << "\nEntering stage 1\n";
        cout << "================\n" << endl;

        if (!settings.stageOneNecessary()) {
                cout << "Files " << settings.getStage1GraphFilename()
                     << ".{node,edge,meta}.st1 exist. Skipping stage 1...\n";
                return;
        }

        Util::startChrono();

        // get all k-mers
        Util::startChrono();
        KmerTable readParser(settings);
        cout << "Generating k-mers (k = " << Kmer::getK()
             << ") from input files..." << endl;
        readParser.parseInputFiles(libraries);
        cout << "Input files contain "
             << readParser.getNumKmers() << " k-mers\n";
        readParser.writeAllKmers("kmers.st1");
        readParser.clear();
        cout << "Generating k-mers in " << Util::stopChronoStr() << "\n";

        // create a kmer table from the reads
        Util::startChrono();
        KmerOverlapTable overlapTable(settings);
        cout << "Building compacted de Bruijn graph..." << endl;
        overlapTable.loadKmersFromDisc("kmers.st1");
        overlapTable.parseInputFiles(libraries);
        overlapTable.extractNodes(settings.getStage1GraphFilename());
        overlapTable.clear();
        cout << "Constructed compacted de Bruijn graph in "
             << Util::stopChronoStr() << "\n";

        cout << "Stage 1 finished in " << Util::stopChronoStr() << endl;
}

void stageTwo(Settings& settings, LibraryContainer& libraries)
{
        cout << "\nEntering stage 2\n";
        cout << "================\n" << endl;

        if (!settings.stageTwoNecessary()) {
                cout << "Files " << settings.getStage2NodeModelFilename()
                     << " and " << settings.getStage2EdgeModelFilename()
                     << " exist. Skipping stage 2...\n";
		return;
        }

        Util::startChrono();
        DBGraph dBG(settings);

        // load the graph from file
        dBG.loadGraph(settings.getStage1GraphFilename());

        // create a <kmer, nodePosPair> table
        cout << "Build a k-mer table..." << endl;
        KmerNPPTable table(dBG, true);
        dBG.createKmerNPPTable(table);

        // stream the reads to get the coverage
        cout << "Compute coverage for nodes and arcs..." << endl;
        dBG.getCovFromReads(libraries, table);
        dBG.writeBinary(settings.getStage2GraphFilename());

        // TODO: uniform initialization of weights is best ??
        vector<double> wN(settings.getNumComp(), 1.0);
        vector<double> wE(settings.getNumComp(), 1.0);

        double initErrCov = settings.getInitErrCov();
        double initCov = settings.getInitCov();
        double initODF = settings.getInitODF();

        if (initCov < 0.0)
                initCov = dBG.getInitialKmerCovEstimate(initErrCov, 0.01);
        cout << "Initial coverage estimate: " << fixed << initCov << endl;

        CovModel edgeModel(initErrCov, initODF, initCov, initODF, wE);
        CovModel nodeModel(initErrCov, initODF, initCov, initODF, wN);

        // train and write coverage models
        trainModel(dBG, settings, nodeModel, edgeModel, true);
        nodeModel.write(settings.getStage2NodeModelFilename());
        edgeModel.write(settings.getStage2EdgeModelFilename());

        cout << "Stage 2 finished in " << Util::stopChronoStr() << endl;
}

NodeMap<int> trueNodeMult; EdgeMap<int> trueEdgeMult;   // FIXME

void stageThreeCorrectGraph(Settings& settings)
{
        cout << "\nEntering stage 3\n";
        cout << "================\n" << endl;

        if (Util::fileExists(settings.getStage3GraphFilename())) {
                cout << "File " << settings.getStage3GraphFilename()
                     << " exist. Skipping stage 3...\n";
		return;
        }

        Util::startChrono();

        // load stage 1 graph
        DBGraph dBG(settings);
        string ifn = settings.getStage2GraphFilename();
        cout << "Loading graph from file: " << ifn << "..."; cout.flush();
        Util::startChrono();
        dBG.loadBinary(ifn);
        cout << "\n\tLoaded " << dBG.getNumNodes() << " nodes and "
             << dBG.getNumArcs() << " arcs (" << Util::stopChronoStr() << ")\n";

        // load node/edge coverage model
        CovModel nodeModel(settings.getStage2NodeModelFilename());
        CovModel edgeModel(settings.getStage2EdgeModelFilename());

        cout << "Node model: " << nodeModel << endl;
        cout << "Edge model: " << edgeModel << endl;

        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMargin(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.getNumThreads(),
                          settings.getCRFWorkSize());

        cout.precision(3);

        /*KmerNPPTable table(dBG, true);
        dBG.createKmerNPPTable(table);

        RefComp refComp(dBG, settings, table, "previous.fasta");
        vector<vector<AlnSegment>> aln;
        refComp.alignSequences(aln);
        refComp.getTrueMultiplicity(aln, trueNodeMult, trueEdgeMult);

        for (auto [nr, m] : trueNodeMult) {
                if (m == 0)
                        continue;
                SSNode n = dBG.getSSNode(nr);
                //if (n.numLeftArcs() == 0 || n.numRightArcs() == 0)
                 //       continue;       // skip tips

                double cov = 10.683 * m * n.getMarginalLength();
                cov = max<double>(cov, n.getCov());
                n.setCov(cov);
        }
        for (auto [er, m] : trueEdgeMult) {
                if (m == 0)
                        continue;
                double cov = 10.975 * m;
                cov = max<double>(cov, dBG.getArc(er).getCov());
                SSNode n = dBG.getSSNode(er.getSrcID());
                //if (n.numLeftArcs() == 0 || n.numRightArcs() == 0)
                //        continue;       // skip tips
                n = dBG.getSSNode(er.getDstID());
                //if (n.numLeftArcs() == 0 || n.numRightArcs() == 0)
                //        continue;       // skip tips

                dBG.setArcCov(er.getSrcID(), er.getDstID(), cov);
        }*/

        // ================
        /*if (Util::fileExists("genome.fasta")) {
                // create a <kmer, nodePosPair> table
                KmerNPPTable table(dBG, true);
                dBG.createKmerNPPTable(table);

                RefComp refComp(dBG, settings, table, "genome.fasta");
                refComp.writeAlignedSeqs("theoretical_BC.fasta");
                refComp.getTrueMultiplicity(trueNodeMult, trueEdgeMult);

                NodeChain nc = refComp.getNodeChain(0, 1701600, 1701650);
                //NodeChain nc = refComp.getNodeChain(0, 2783320, 2783430);
                refComp.getTrueMultiplicity(trueNodeMult, trueEdgeMult);
                for (NodeID id : nc)
                        cout << id << (trueNodeMult[id] == 1 ? " (U)\t" : "\t");
                cout << endl;

                vector<NodeID> nodes; vector<EdgeID> edges;
                if (!nc.empty())
                        dBG.getSubgraph(nc, nodes, edges, 1);
                //dBG.getSubgraph(2122756, nodes, edges, 3);

                cout << "Wrote Cytoscape graph (" << nodes.size() << " nodes)" << endl;

                NodeMap<Multiplicity> nodeMult(nodes.size());
                EdgeMap<Multiplicity> edgeMult(edges.size());
                for (size_t i = 0; i < nodes.size(); i++)
                        nodeMult[nodes[i]] = Multiplicity();
                for (size_t i = 0; i < edges.size(); i++)
                        edgeMult[edges[i]] = Multiplicity();
                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel);

                dBG.writeCytoscapeGraph("graph", nodes, edges, nodeMult, edgeMult,
                                        trueNodeMult, trueEdgeMult);
                //exit(0);
        }*/
        // ================

        double maxCov = nodeModel.getCovCutOff(1);
        vector<double> covThreshold = nodeModel.getFractions(maxCov);

        for (size_t i = 0; i < covThreshold.size() + 1; i++)
        {
                double maxCov = (i < covThreshold.size()) ?
                        covThreshold[i] : covThreshold.back();

                size_t numRmvd = dBG.removeLowCovTips(maxCov, MAX_TIP_LEN);
                cout << "Removed " << numRmvd << " tips "
                     << "with coverage <= " << maxCov << "\n";
                cout << "\tGraph has " << dBG.getNumValidNodes()
                     << " nodes and " << dBG.getNumValidArcs() << " arcs\n";

                numRmvd = dBG.removeLowCovBubbles(maxCov, MAX_BUBBLE_LEN);
                cout << "Removed " << numRmvd << " bubbles "
                     << "with coverage <= " << maxCov << "\n";
                cout << "\tGraph has " << dBG.getNumValidNodes()
                     << " nodes and " << dBG.getNumValidArcs() << " arcs\n";

                dBG.removeLowCovNodesFlow(maxCov, myCRFMult,
                                          nodeModel, edgeModel);
                dBG.removeLowCovEdgesFlow(maxCov, myCRFMult,
                                          nodeModel, edgeModel);

                dBG.removeLowCovNodesCRF(maxCov, myCRFMult,
                                         nodeModel, edgeModel);
                dBG.removeLowCovEdgesCRF(maxCov, myCRFMult,
                                         nodeModel, edgeModel);

                if (i+1 < covThreshold.size())
                        continue;

                // Remove high-coverage tips if dictated by the CRF
                do {
                        numRmvd = dBG.removeTipsCRF(MAX_TIP_LEN, myCRFMult, nodeModel, edgeModel);
                        cout << "Removed " << numRmvd << " high-coverage tips\n";
                } while (numRmvd > 0);

                // Remove high-coverage bubbles if dictated by the CRF
                do {
                        numRmvd = dBG.removeBubblesCRF(MAX_TIP_LEN, myCRFMult, nodeModel, edgeModel);
                        cout << "Removed " << numRmvd << " high-coverage bubbles\n";
                } while (numRmvd > 0);

                //exit(0);

                // high coverage bubbles
                /*dBG.removeLowCovBubbles(1e100, MAX_BUBBLE_LEN);
                cout << "Removed " << numRmvd
                     << " bubbles with coverage <= " << 1e100 << endl;
                dBG.concatenateNodes();
                cout << "\tGraph has " << dBG.getNumValidNodes()
                     << " nodes and " << dBG.getNumValidArcs() << " arcs\n";*/

                // High coverage tips
                /*vector <NodeRep> nodeReps = dBG.getLowCovTips(1e100);
                cout << "Selected " << nodeReps.size() << " tips" << endl;
                NodeMap<Multiplicity> nodeMult(nodeReps.size());
                for (size_t i = 0; i < nodeReps.size(); i++)
                        nodeMult[nodeReps[i]] = Multiplicity();
                EdgeMap<Multiplicity> edgeMult;
                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel);
                removeNodes(dBG, nodeMult);*/
        }

        // ================ DEBUG ================
        /*NodeMap<int> trueNodeMult; EdgeMap<int> trueEdgeMult;
        if (Util::fileExists("genome.fasta")) {
                // create a <kmer, nodePosPair> table
                KmerNPPTable table(dBG, true);
                dBG.createKmerNPPTable(table);

                RefComp refComp(dBG, settings, table, "genome.fasta");
                refComp.writeAlignedSeqs("theoretical_AC.fasta");
                refComp.getTrueMultiplicity(trueNodeMult, trueEdgeMult);

                //NodeChain nc = refComp.getNodeChain(0, 58180, 58290);
                NodeChain nc = refComp.getNodeChain(0, 1701600, 1701650);
                //NodeChain nc = refComp.getNodeChain(0, 2783320, 2783430);
                refComp.getTrueMultiplicity(trueNodeMult, trueEdgeMult);
                for (NodeID id : nc)
                        cout << id << (trueNodeMult[id] == 1 ? " (U)\t" : "\t");
                cout << endl;

                vector<NodeID> nodes; vector<EdgeID> edges;
                if (!nc.empty())
                        dBG.getSubgraph(nc, nodes, edges, 2);
                //dBG.getSubgraph(2122756, nodes, edges, 3);

                cout << "Wrote Cytoscape graph (" << nodes.size() << " nodes)" << endl;

                NodeMap<Multiplicity> nodeMult(nodes.size());
                EdgeMap<Multiplicity> edgeMult(edges.size());
                for (size_t i = 0; i < nodes.size(); i++)
                        nodeMult[nodes[i]] = Multiplicity();
                for (size_t i = 0; i < edges.size(); i++)
                        edgeMult[edges[i]] = Multiplicity();
                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel);
                dBG.writeCytoscapeGraph("graph", nodes, edges, nodeMult, edgeMult,
                                        trueNodeMult, trueEdgeMult);
                //exit(0);
        }*/
        // ================ DEBUG ================

        // re-train and write coverage models (do NOT retrain node ODF)
        trainModel(dBG, settings, nodeModel, edgeModel, false);
        nodeModel.write(settings.getStage3NodeModelFilename());
        edgeModel.write(settings.getStage3EdgeModelFilename());

        dBG.defrag();
        dBG.sanityCheck();

        // write stage 3 binary output file
        dBG.writeBinary(settings.getStage3GraphFilename());
        cout << "Stage 3 finished in " << Util::stopChronoStr() << endl;
}

void stageFourResolveRepeats(Settings& settings, LibraryContainer& libraries)
{
        cout << "\nEntering stage 4\n";
        cout << "================\n" << endl;\

        if (Util::fileExists("detox.st4.fasta")) {
                cout << "File " << "detox.st4.fasta"
                     << " exist. Skipping stage 4...\n";
		return;
        }

        Util::startChrono();

        // load stage 3 graph
        string ifn = settings.getStage3GraphFilename();
        DBGraph dBG(settings);
        cout << "Loading graph from file: " << ifn << "..."; cout.flush();
        Util::startChrono();
        dBG.loadBinary(ifn);
        cout << "\n\tLoaded " << dBG.getNumNodes() << " nodes and "
             << dBG.getNumArcs() << " arcs (" << Util::stopChronoStr() << ")\n";

        // reset all flags
        dBG.setAllFlags1(false);
        dBG.setAllFlags2(false);

        // load the CRF and, optionally, neural network models
        MultEstimation multEst(dBG, settings,
                               settings.getStage3NodeModelFilename(),
                               settings.getStage3EdgeModelFilename(),
                               settings.getNNModelFilename());

        // align reads to graph
        KmerNPPTable table(dBG, true);
        dBG.createKmerNPPTable(table);
        ReadAlnHandler raHandler(dBG, settings, table);
        Util::startChrono();
        raHandler.align(libraries);
        libraries.filterAlignments(dBG);
        cout << "Aligned reads in " << Util::stopChronoStr() << endl;

        // get the true node multiplicities
        /*RefComp refComp(dBG, settings, table, "genome.fasta");
        refComp.writeAlignedSeqs("theoretical_AC.fasta");
        NodeMap<int> trueNodeMult; EdgeMap<int> trueEdgeMult;
        refComp.getTrueMultiplicity(trueNodeMult, trueEdgeMult);*/

        for (int i = 1; i <= 1; i++) {
                cout << "=============== PASS " << i << " ===============\n";

                // A) Read-based reductions

                // estimate node/edge multiplicities using the CRF model
                auto [nodeMult, edgeMult] = multEst.computeMultiplicities();

                // using the estimated multiplicities, get the reductions
                vector<NodeChain> reductions, looseEnds;
                dBG.getReductions(libraries, nodeMult, edgeMult,
                                  settings.getMCD(), settings.getMinCount(),
                                  reductions, looseEnds);

                // adjust estimated multiplicity according to the reductions
                NodeMap<int> redNodeMult; EdgeMap<int> redEdgeMult;
                dBG.adjustMult(reductions, looseEnds, nodeMult, edgeMult,
                               redNodeMult, redEdgeMult);

                // build the node tracker
                NodeTracker nodeTracker(dBG);
                nodeTracker.build(reductions, nodeMult, redNodeMult);

                // apply the reductions in the de Bruijn graph
                cout << "Applying " << reductions.size() << " read based reductions\n";
                dBG.applyReductions(reductions, nodeMult, edgeMult);
                cout << "\tGraph has " << dBG.getNumValidNodes() << " nodes and "
                     << dBG.getNumValidArcs() << " arcs" << endl;

                // adapt the read alignments according to the reductions
                libraries.parseAlignments(nodeTracker);

                // ============================================================

                // B) Concatenations

                // get the concatenations
                dBG.getConcatenations(reductions);

                // make sure the concatenated node get removed
                for (const auto& r : reductions) {
                        for (size_t i = 0; i < r.size(); i++)
                                nodeMult[NodeRep(r[i])] = 1;
                        for (size_t i = 1; i < r.size(); i++)
                                edgeMult[EdgeRep(r[i-1], r[i])] = 1;
                }

                // apply the reductions in the de Bruijn graph
                cout << "Applying " << reductions.size() << " concatenations\n";
                dBG.applyReductions(reductions, nodeMult, edgeMult);
                cout << "\tGraph has " << dBG.getNumValidNodes() << " nodes and "
                     << dBG.getNumValidArcs() << " arcs" << endl;

                // adapt the read alignments according to the reductions
                nodeTracker.build(reductions);
                libraries.parseAlignments(nodeTracker);

                // ============================================================

                // C) Graph-based reductions
                tie(nodeMult, edgeMult) = multEst.computeMultiplicities();

                vector<NodePair> tangles;
                dBG.findTangles(nodeMult, tangles);

                dBG.resolveTangles(libraries, nodeMult, edgeMult, tangles,
                                   reductions, redNodeMult, redEdgeMult);

                cout << "Applying " << reductions.size() << " graph-based reductions\n";
                dBG.applyReductions(reductions, redNodeMult, redEdgeMult);
                cout << "\tGraph has " << dBG.getNumValidNodes() << " nodes and "
                     << dBG.getNumValidArcs() << " arcs" << endl;

                // adapt the read alignments according to the reductions
                nodeTracker.build(reductions);
                libraries.parseAlignments(nodeTracker);

                // ============================================================

                // D) Linear unique paths
                //computeAllMultiplicities(dBG, settings, nodeMult, edgeMult);

                //dBG.findLinearUniquePaths(nodeMult, edgeMult);
        }

        // =================================================
        dBG.glueTips(libraries);
        // =================================================

        // get the node multiplicities
        CovModel nodeModel(settings.getStage3NodeModelFilename());
        CovModel edgeModel(settings.getStage3EdgeModelFilename());

        cout << "Node model: " << nodeModel << endl;
        cout << "Edge model: " << edgeModel << endl;

        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMargin(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.getNumThreads(),
                          settings.getCRFWorkSize());

        // Remove high-coverage bubbles if dictated by the CRF

        while (dBG.removeTipsCRF(MAX_TIP_LEN, myCRFMult, nodeModel, edgeModel) > 0);
        while (dBG.removeBubblesCRF(MAX_TIP_LEN, myCRFMult, nodeModel, edgeModel) > 0);
        // =================================================

        NodeMap<Multiplicity> nodeMult; EdgeMap<Multiplicity> edgeMult;
        computeAllMultiplicities(dBG, settings, nodeMult, edgeMult);

        // write final set of contigs
        vector<NodeChain> contigs; vector<string> contigNames;
        //dBG.getContigs(libraries, nodeMulti, edgeMulti, contigs, contigNames);
        dBG.getSingletonContigs(nodeMult, contigs, contigNames);

        dBG.writeContigs("detox.fasta", contigs, contigNames);

        // write Cytoscape graph
        vector<NodeID> nodes; vector<EdgeID> edges;
        /*nc = vector<NodeID>({9154, -8497, 1069});
        if (!nc.empty())
               dBG.getSubgraph(nc, nodes, edges, 4);*/
        dBG.getGraph(nodes, edges);
        //dBG.getSubgraph(9154, nodes, edges, 5);
        cout << "Wrote Cytoscape graph (" << nodes.size() << " nodes)" << endl;
        dBG.writeCytoscapeGraph("graph", nodes, edges, nodeMult, edgeMult,
                                trueNodeMult, trueEdgeMult);

        // write stage 4 binary output file
        dBG.defrag();
        dBG.writeBinary(settings.getStage4GraphFilename());

        cout << "Stage 4 finished in " << Util::stopChronoStr() << endl;
}

int main(int argc, char** argv)
{
        try {
                Settings settings(argc, argv);
                LibraryContainer libraries(settings.getReadFilename());

                cout << "Welcome to Detox version "
                     << DETOX_MAJOR_VERSION << "."
                     << DETOX_MINOR_VERSION << "."
                     << DETOX_PATCH_LEVEL;

#ifdef DEBUG
                cout << " (debug mode)" << endl;
#else
                cout << " (release mode)" << endl;
#endif
                stageOne(settings, libraries);
                stageTwo(settings, libraries);
                stageThreeCorrectGraph(settings);
                stageFourResolveRepeats(settings, libraries);
        } catch (exception& e) {
                cerr << "Fatal error: " << e.what() << endl;
                return EXIT_FAILURE;
        }

        cout << "Exiting... bye!" << endl;
        return EXIT_SUCCESS;
}
