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

#include <cstdlib>
#include <numeric>

#include "dbgraph.h"
#include "library.h"
#include "settings.h"
#include "refcomp.h"
#include "crfmult.h"

using namespace std;

void stageOne(Settings& settings, LibraryContainer& libraries)
{
        cout << "\nEntering stage 1\n";
        cout << "================\n" << endl;

        if (!settings.stageOneNecessary()) {
                cout << "File " << settings.getStage1GraphFilename()
                     << " exists. Skipping stage 1...\n";
                return;
        }

        Util::startChrono();

        DBGraph dBG(settings);

        // load the graph from BCALM2 file
        dBG.loadBCalm(settings.getGraphFilename());

        // create a <kmer, nodePosPair> table
        KmerNPPTable table(dBG, true);
        dBG.createKmerNPPTable(table);

        // stream the reads to get the coverage
        dBG.getCovFromReads(libraries, table);

        // write stage 1 binary output file
        dBG.writeBinary(settings.getStage1GraphFilename());

        if (Util::fileExists("genome.fasta")) {
                //dBG.sanityCheck();

                RefComp refComp(dBG, settings, "genome.fasta");
                //refComp.alignSequences();

                vector<int> trueNodeMult, trueEdgeMult;
                refComp.getTrueMultiplicity(table, trueNodeMult, trueEdgeMult);

                vector<NodeRep> nodes = dBG.getNodeReps(dBG.getNumValidNodes());
                vector<EdgeRep> edges = dBG.getEdgeReps(dBG.getNumValidArcs());

                ofstream nodeOFS(settings.getTrueNodeMultFilename());
                ofstream edgeOFS(settings.getTrueEdgeMultFilename());

                for (NodeID id = 1; id <= dBG.getNumNodes(); id++) {
                        SSNode n = dBG.getSSNode(id);
                        if (!n.isValid())
                                continue;

                        nodeOFS << id << "\t" << trueNodeMult[id] << "\n";

                        for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++)
                                if (id <= -it->getNodeID()) {
                                        EdgeRep eRep(id, it->getNodeID());
                                        edgeOFS << eRep.getSrcID() << "\t" << eRep.getDstID()
                                                << "\t" << trueEdgeMult[dBG.getArcID(eRep)] << "\n";
                                }
                        for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++)
                                if (-id <= it->getNodeID()) {
                                        EdgeRep eRep(-id, -it->getNodeID());
                                        edgeOFS << eRep.getSrcID() << "\t" << eRep.getDstID()
                                                << "\t" << trueEdgeMult[dBG.getArcID(eRep)] << "\n";

                                }
                }

                nodeOFS.close();
                edgeOFS.close();
        }

        cout << "Stage 1 finished in " << Util::stopChronoStr() << endl;
}

void stageTwo(Settings& settings)
{
        cout << "\nEntering stage 2\n";
        cout << "================\n" << endl;

        if (!settings.stageTwoNecessary()) {
                cout << "Files " << settings.getNodeModelFilename() << " and "
                     << settings.getEdgeModelFilename() << " exist. Skipping stage 2...\n";
		return;
        }

        Util::startChrono();

        // load stage 1 de Bruijn graph
        DBGraph dBG(settings);
        dBG.loadBinary(settings.getStage1GraphFilename());

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

        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMargin(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.getNumThreads(),
                          settings.getThreadGraphWorkSize());

        vector<NodeRep> nodes = dBG.getNodeReps(settings.getEMTrainingSize());
        vector<EdgeRep> edges = dBG.getEdgeReps(settings.getEMTrainingSize());
        vector<Multiplicity> nodeMult(nodes.size()), edgeMult(edges.size());

        cout.precision(2);
        cout << "Fitting mixture model using EM for " << nodes.size()
             << " nodes and " << edges.size() << " edges." << endl;

        int nIters = myCRFMult.computeMultEM(nodes, nodeMult, nodeModel,
                                             edges, edgeMult, edgeModel,
                                             settings.getEMConvEps(),
                                             settings.getEMMaxIter());

        if (nIters <= settings.getEMMaxIter())
                cout << "EM algorithm converged after " << nIters << " iterations\n";
        else
                cout << "WARNING: maximum number of iterations reached. "
                     << "Convergence is not guaranteed." << endl;

        // write fitted coverage models
        nodeModel.write(settings.getNodeModelFilename());
        edgeModel.write(settings.getEdgeModelFilename());

        map<int, double> nodeHist;
        for (size_t i = 0; i < nodes.size(); i++) {
                SSNode node = dBG.getSSNode(nodes[i]);
                double f = node.getAvgCov() - floor(node.getAvgCov());
                nodeHist[node.getAvgCov()] += (1.0-f);
                nodeHist[node.getAvgCov() + 1] += f;
        }

        nodeModel.writeGnuplot("nodes", nodeHist);

        map<int, double> edgeHist;
        for (size_t i = 0; i < edges.size(); i++) {
                Arc& arc = dBG.getArc(dBG.getArcID(edges[i]));
                double f = arc.getCov() - floor(arc.getCov());
                edgeHist[arc.getCov()] += (1.0-f);
                edgeHist[arc.getCov() + 1] += f;
        }

        edgeModel.writeGnuplot("edges", edgeHist);

        cout << "Stage 2 finished in " << Util::stopChronoStr() << endl;
}

void stageThreeCompAllMult(Settings& settings, DBGraph& dBG, CovModel& nodeModel,
                           CovModel& edgeModel, CRFMult& myCRFMult)
{
        Util::startChrono();

        cout << "Inferring node multiplicities for " << dBG.getNumValidNodes()
             << " nodes and " << dBG.getNumValidArcs() << " edges...\n";
        cout << "Using CRF neighborhood size: " << settings.getCRFDepth()
             << " and using " << settings.getNumThreads() << " threads... " << endl;;

        // output files (estimated multiplicities for nodes/edges)
        ofstream nodeOFS(settings.getEstNodeMultFilename());
        ofstream edgeOFS(settings.getEstEdgeMultFilename());

        WorkLoadBalancer wlb(1, dBG.getNumNodes()+1, 10000000, "");
        size_t idxBegin, idxEnd;
        while (wlb.getChunk(idxBegin, idxEnd)) {
                cout << "\tHandling nodes " << idxBegin << " to " << idxEnd-1 << "...\n";

                // prepare the next chunk of nodes and edges
                vector<NodeRep> nodes; vector<EdgeRep> edges;
                nodes.reserve(idxEnd - idxBegin);
                edges.reserve(idxEnd - idxBegin);

                for (NodeID id = idxBegin; id < idxEnd; id++) {
                        SSNode n = dBG.getSSNode(id);
                        if (n.isValid())
                                nodes.push_back(NodeRep(id));
                        for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++)
                                if (id <= -it->getNodeID())
                                        edges.push_back(EdgeRep(id, it->getNodeID()));
                        for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++)
                                if (-id <= it->getNodeID())
                                        edges.push_back(EdgeRep(-id, -it->getNodeID()));
                }

                vector<Multiplicity> nodeMult(nodes.size()), edgeMult(edges.size());
                myCRFMult.computeMult(nodes, nodeMult, nodeModel,
                                      edges, edgeMult, edgeModel);

                for (size_t i = 0; i < nodeMult.size(); i++) {
                        SSNode node = dBG.getSSNode(nodes[i]);
                        nodeOFS << nodes[i].getNodeID() << "\t"
                                << nodeMult[i].getExpMult() << "\t"
                                << nodeMult[i].getExpMultLProb() << "\t"
                                << node.getAvgCov() << "\t"
                                << node.getMarginalLength() << "\n";
                }

                for (size_t i = 0; i < edgeMult.size(); i++) {
                        edgeOFS << edges[i].getSrcID() << "\t"
                                << edges[i].getDstID() << "\t"
                                << edgeMult[i].getExpMult() << "\t"
                                << edgeMult[i].getExpMultLProb() << "\t"
                                << dBG.getArc(edges[i]).getCov() << "\n";
                }
        }

        cout << "Wrote file: " << settings.getEstNodeMultFilename() << endl;
        cout << "Wrote file: " << settings.getEstEdgeMultFilename() << endl;

        CRFPerfCounter pc = myCRFMult.getPerfCounter();
        cout << "Performance metrics:\n";
        cout << "\tAvg. no. of nodes in the CRF neighborhood: " << pc.getAvgNumNodes() << endl;
        cout << "\tAvg. no. of edges in the CRF neighborhood: " << pc.getAvgNumEdges() << endl;
        cout << "\tMultiplicity computation finished in: " << Util::stopChronoStr() << endl;
}

void stageThreeCompMultFile(Settings& settings, DBGraph& dBG, CovModel& nodeModel,
                           CovModel& edgeModel, CRFMult& myCRFMult)
{
        Util::startChrono();

        cout << "Inferring node multiplicities for nodes and edges stored in files "
             << settings.getNodeListFilename() << " and " << settings.getEdgeListFilename() << "\n";
        cout << "Using CRF neighborhood size: " << settings.getCRFDepth()
             << " and using " << settings.getNumThreads() << " threads... " << endl;;

        // input files (list of nodes/edges for which to compute the multiplicity)
        ifstream nodeIFS(settings.getNodeListFilename());
        ifstream edgeIFS(settings.getEdgeListFilename());

        // output files (estimated multiplicities for nodes/edges)
        ofstream nodeOFS(settings.getEstNodeMultFilename());
        ofstream edgeOFS(settings.getEstEdgeMultFilename());

        size_t chunkSize = 10000000;
        while (nodeIFS || edgeIFS) {
                cout << "\tHandling next chunk of nodes and edges...\n";

                // prepare the next chunk of nodes and edges
                vector<NodeRep> nodes; vector<EdgeRep> edges;
                nodes.reserve(chunkSize);
                edges.reserve(chunkSize);

                for (size_t i = 0; i < chunkSize; i++) {
                        NodeID nodeID;
                        if (nodeIFS >> nodeID)
                                if (dBG.nodeExists(nodeID))
                                        nodes.push_back(nodeID);

                        ArcID srcID, dstID;
                        if (edgeIFS >> srcID >> dstID)
                                if (dBG.edgeExists(EdgeRep(srcID, dstID)))
                                        edges.push_back(EdgeRep(srcID, dstID));
                }

                vector<Multiplicity> nodeMult(nodes.size()), edgeMult(edges.size());
                myCRFMult.computeMult(nodes, nodeMult, nodeModel,
                                      edges, edgeMult, edgeModel);

                for (size_t i = 0; i < nodeMult.size(); i++) {
                        SSNode node = dBG.getSSNode(nodes[i]);
                        nodeOFS << nodes[i].getNodeID() << "\t"
                                << nodeMult[i].getExpMult() << "\t"
                                << nodeMult[i].getExpMultLProb() << "\t"
                                << node.getAvgCov() << "\t"
                                << node.getMarginalLength() << "\n";
                }

                for (size_t i = 0; i < edgeMult.size(); i++) {
                        edgeOFS << edges[i].getSrcID() << "\t"
                                << edges[i].getDstID() << "\t"
                                << edgeMult[i].getExpMult() << "\t"
                                << edgeMult[i].getExpMultLProb() << "\t"
                                << dBG.getArc(edges[i]).getCov() << "\n";
                }
        }

        cout << "Wrote file: " << settings.getEstNodeMultFilename() << endl;
        cout << "Wrote file: " << settings.getEstEdgeMultFilename() << endl;

        CRFPerfCounter pc = myCRFMult.getPerfCounter();
        cout << "Performance metrics:\n";
        cout << "\tAvg. no. of nodes in the CRF neighborhood: " << pc.getAvgNumNodes() << endl;
        cout << "\tAvg. no. of edges in the CRF neighborhood: " << pc.getAvgNumEdges() << endl;
        cout << "\tMultiplicity computation finished in: " << Util::stopChronoStr() << endl;
}

void stageThreeCompMult(Settings& settings)
{
        cout << "\nEntering stage 3\n";
        cout << "================\n" << endl;

        Util::startChrono();

        DBGraph dBG(settings);
        dBG.loadBinary(settings.getStage1GraphFilename());

        cout << "Loaded graph with " << dBG.getNumValidNodes() << " nodes and "
             << dBG.getNumValidArcs() << " edges." << endl;

        CovModel nodeModel(settings.getNodeModelFilename());
        CovModel edgeModel(settings.getEdgeModelFilename());

        cout << "Node model: " << nodeModel << endl;
        cout << "Edge model: " << edgeModel << endl;

        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMargin(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.getNumThreads(),
                          settings.getThreadGraphWorkSize());

        if ( Util::fileExists(settings.getNodeListFilename()) ||
             Util::fileExists(settings.getEdgeListFilename()) )
                stageThreeCompMultFile(settings, dBG, nodeModel, edgeModel, myCRFMult);
        else
                stageThreeCompAllMult(settings, dBG, nodeModel, edgeModel, myCRFMult);

        cout << "Stage 3 finished in " << Util::stopChronoStr() << endl;
}

void stageThreeCG(Settings& settings)
{
        cout << "\nEntering stage 3\n";
        cout << "================\n" << endl;

        DBGraph dBG(settings);
        dBG.loadBinary(settings.getStage1GraphFilename());

        NodeID noi = settings.getCytGraphNode();
        if ((noi < -dBG.getNumNodes()) || (noi > dBG.getNumNodes()))
                throw runtime_error("Specified node id of " + to_string(noi) +
                                    " does not exist in the graph");

        CovModel nodeModel(settings.getNodeModelFilename());
        CovModel edgeModel(settings.getEdgeModelFilename());

        cout << "Node model: " << nodeModel << endl;
        cout << "Edge model: " << edgeModel << endl;

        vector<NodeID> nodes;
        vector<pair<NodeID, NodeID> > edges;
        dBG.getSubgraph(noi, nodes, edges, settings.getCRFDepth());

        vector<NodeRep> nodeReps(nodes.begin(), nodes.end());
        vector<EdgeRep> edgeReps(edges.begin(), edges.end());

        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMargin(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.getNumThreads(),
                          settings.getThreadGraphWorkSize());

        vector<Multiplicity> nodeMult(nodes.size()), edgeMult(edges.size());

        myCRFMult.computeMult(nodeReps, nodeMult, nodeModel,
                              edgeReps, edgeMult, edgeModel);

        // Convert estimated node/edge multiplicities to lookup table
        map<NodeRep, Multiplicity> estNodeMult;
        map<EdgeRep, Multiplicity> estEdgeMult;

        for (size_t i = 0; i < nodeReps.size(); i++)
                estNodeMult[nodeReps[i]] = nodeMult[i];
        for (size_t i = 0; i < edgeReps.size(); i++)
                estEdgeMult[edgeReps[i]] = edgeMult[i];

        // If they exist, also load the true multiplicities from disk
        // Note: for efficiency, story only true multiplicities that are needed
        map<NodeRep, int> trueNodeMult;
        transform(nodeReps.begin(), nodeReps.end(),
                  inserter(trueNodeMult, trueNodeMult.end()),
                  [](const NodeRep& nr) { return make_pair(nr, Multiplicity(-1)); });

        map<EdgeRep, int> trueEdgeMult;
        transform(edgeReps.begin(), edgeReps.end(),
                  inserter(trueEdgeMult, trueEdgeMult.end()),
                  [](const EdgeRep& er) { return make_pair(er, Multiplicity(-1)); });

        ifstream ifs("truemult.node.st1");
        if (!ifs)
                cerr << "Could not find true nodes multiplicities file...\n"
                        "True node multiplicities will be set to -1 in "
                        "resulting Cytoscape graph\n";

        while (ifs) {
                NodeID nodeID; int m;
                ifs >> nodeID >> m;
                if (!ifs)
                        break;
                auto it = trueNodeMult.find(nodeID);
                if (it != trueNodeMult.end())
                        it->second = m;
        }
        ifs.close();

        ifs.open("truemult.arc.st1");
        if (!ifs)
                cerr << "Could not find true edges multiplicities file...\n"
                        "True edge multiplicities will be set to -1 in "
                        "resulting Cytoscape graph\n";

        while (ifs) {
                NodeID srcID, dstID; int m;
                ifs >> srcID >> dstID >> m;
                if (!ifs)
                        break;
                auto it = trueEdgeMult.find(EdgeRep(srcID, dstID));
                if (it != trueEdgeMult.end())
                        it->second = m;
        }
        ifs.close();

        cout << "Writing Cytoscape graph... " << endl;
        dBG.writeCytoscapeGraph("cytgraph" + to_string(noi), nodes, edges,
                                estNodeMult, estEdgeMult,
                                trueNodeMult, trueEdgeMult);
}

void stageThreeCorrectGraph(Settings& settings)
{
        cout << "\nEntering stage 3\n";
        cout << "================\n" << endl;

        Util::startChrono();

        DBGraph dBG(settings);

        // load stage 1 graph
        string ifn = settings.getStage1GraphFilename();
        cout << "Loading graph from file: " << ifn << "..."; cout.flush();
        Util::startChrono();
        dBG.loadBinary(settings.getStage1GraphFilename());
        cout << "\n\tLoaded " << dBG.getNumNodes() << " nodes and "
             << dBG.getNumArcs() << " arcs (" << Util::stopChronoStr() << ")\n";

#ifdef DEBUG
        {
        RefComp refComp(dBG, settings, "genome.fasta");
        refComp.buildKmerNPPTable();
        refComp.alignSequences();

        vector<int> trueNodeMult, trueEdgeMult;
        refComp.getTrueMultiplicity(trueNodeMult, trueEdgeMult);
        }
#endif

        // load node/edge coverage model
        CovModel nodeModel(settings.getNodeModelFilename());
        CovModel edgeModel(settings.getEdgeModelFilename());

        cout << "Node model: " << nodeModel << endl;
        cout << "Edge model: " << edgeModel << endl;

        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMargin(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.getNumThreads(),
                          settings.getThreadGraphWorkSize());

        // remove arcs with zero coverage
        cout << "Removing nodes/arcs with zero coverage..." << endl;
        dBG.removeCoverage(0.0, 0);
        dBG.concatenateNodes();

        cout << "\tGraph has " << dBG.getNumValidNodes() << " nodes and "
             << dBG.getNumValidArcs() << " arcs" << endl;

        // remove low-coverage tips if the flow is OK
        double nodeCutoff = nodeModel.getCovCutOff(1);

        while (true) {
                vector<NodeRep> nodeReps = dBG.getLowCovTips(nodeCutoff);
                cout << "Selected " << nodeReps.size()
                     << " tips with coverage <= " << nodeCutoff << endl;

                vector<bool> flowOK(nodeReps.size());
                myCRFMult.checkFlow(nodeReps, flowOK, nodeModel, edgeModel);

                size_t numRemove = count(flowOK.begin(), flowOK.end(), true);
                if (numRemove == 0)
                        break;

                vector<NodeRep> toRemove;
                toRemove.reserve(numRemove);

                for (size_t i = 0; i < nodeReps.size(); i++)
                        if (flowOK[i])
                                toRemove.push_back(nodeReps[i]);

                cout << "\tRemoving " << numRemove << " nodes" << endl;
                dBG.removeNodes(toRemove);
                dBG.concatenateNodes();

                cout << "\tGraph has " << dBG.getNumValidNodes() << " nodes and "
                     <<  dBG.getNumValidArcs() << " arcs" << endl;
        }

        // remove low-coverage bubbles if the flow is OK
        {
                vector<NodeRep> nodeReps = dBG.getLowCovBubbles(nodeCutoff);
                cout << "Selected " << nodeReps.size()
                     << " bubbles with coverage <= " << nodeCutoff << endl;

                vector<bool> flowOK(nodeReps.size());
                myCRFMult.checkFlow(nodeReps, flowOK, nodeModel, edgeModel);

                size_t numRemove = count(flowOK.begin(), flowOK.end(), true);
                vector<NodeRep> toRemove;
                toRemove.clear();
                toRemove.reserve(numRemove);

                cout << "\tRemoving " << numRemove << " nodes" << endl;
                for (size_t i = 0; i < nodeReps.size(); i++)
                        if (flowOK[i])
                                toRemove.push_back(nodeReps[i]);

                dBG.removeNodes(toRemove);
                dBG.concatenateNodes();

                cout << "\tGraph has " << dBG.getNumValidNodes() << " nodes and "
                     <<  dBG.getNumValidArcs() << " arcs" << endl;
        }

        // remove low-coverage nodes if the CRF is OK
        {
                vector<NodeRep> nodeReps = dBG.getLowCovNodes(nodeCutoff);
                cout << "Selected " << nodeReps.size()
                     << " nodes with coverage <= " << nodeCutoff << endl;

                vector<Multiplicity> nodeMult(nodeReps.size()), edgeMult;
                myCRFMult.computeMult(nodeReps, nodeMult, nodeModel,
                                      vector<EdgeRep>(), edgeMult, edgeModel);

                vector<NodeRep> toRemove;

                for (size_t i = 0; i < nodeReps.size(); i++)
                        if (nodeMult[i].getExpMult() == 0)
                                toRemove.push_back(nodeReps[i]);

                cout << "\tRemoving " << toRemove.size() << " nodes" << endl;
                dBG.removeNodes(toRemove);
                dBG.concatenateNodes();

                cout << "\tGraph has " << dBG.getNumValidNodes() << " nodes and "
                     <<  dBG.getNumValidArcs() << " arcs" << endl;
        }

        // remove nodes if the CRF is OK
        {
                vector<NodeRep> nodeReps = dBG.getNodeReps(dBG.getNumValidNodes());
                vector<EdgeRep> edgeReps = dBG.getEdgeReps(dBG.getNumValidArcs());
                cout << "Selected all " << nodeReps.size() << " nodes" << endl;

                vector<Multiplicity> nodeMult(nodeReps.size());
                vector<Multiplicity> edgeMult(edgeReps.size());
                myCRFMult.computeMult(nodeReps, nodeMult, nodeModel,
                                      edgeReps, edgeMult, edgeModel);

                vector<NodeRep> nodesToRemove;
                for (size_t i = 0; i < nodeReps.size(); i++)
                        if (nodeMult[i].getExpMult() == 0)
                                nodesToRemove.push_back(nodeReps[i]);

                cout << "\tRemoving " << nodesToRemove.size() << " nodes" << endl;
                dBG.removeNodes(nodesToRemove);

                vector<EdgeRep> edgesToRemove;
                for (size_t i = 0; i < edgeReps.size(); i++)
                        if (edgeMult[i].getExpMult() == 0)
                                edgesToRemove.push_back(edgeReps[i]);

                cout << "\tRemoving " << edgesToRemove.size() << " arcs" << endl;
                dBG.removeEdges(edgesToRemove);

                dBG.concatenateNodes();

                cout << "\tGraph has " << dBG.getNumValidNodes() << " nodes and "
                     <<  dBG.getNumValidArcs() << " arcs" << endl;
        }

        dBG.defrag();
        dBG.sanityCheck();

#ifdef DEBUG
        {
        RefComp refComp(dBG, settings, "genome.fasta");
        refComp.buildKmerNPPTable();
        refComp.alignSequences();

        vector<int> trueNodeMult, trueEdgeMult;
        refComp.getTrueMultiplicity(trueNodeMult, trueEdgeMult);
        }
#endif

        // write stage 3 binary output file
        dBG.writeBinary(settings.getStage3GraphFilename());
        ofstream ofs("trusted.fa");
        vector<NodeRep> nodeReps = dBG.getNodeReps(dBG.getNumValidNodes());
        for (size_t i = 0; i < nodeReps.size(); i++) {
                SSNode node = dBG.getSSNode(nodeReps[i].getNodeID());
                ofs << ">" << nodeReps[i].getNodeID() << "\t"
                    << "COV:" << node.getAvgCov() << "\t"
                    << "LN:" << node.getMarginalLength() << "\n"
                    << node.getSequence()
                    << endl;
        }
        ofs.close();

        double elapsed = Util::stopChrono();
        cout << "Stage 3 finished in " << elapsed << endl;

}

int main(int argc, char** argv)
{
        try {
                Settings settings(argc, argv);
                LibraryContainer libraries(settings.getReadFilename());


                cout << "Welcome to Detox v" << DETOX_MAJOR_VERSION << "."
                     << DETOX_MINOR_VERSION << "." << DETOX_PATCH_LEVEL;

#ifdef DEBUG
                cout << " (debug mode)" << endl;
#else
                cout << " (release mode)" << endl;
#endif

                stageOne(settings, libraries);
                stageTwo(settings);

                if (settings.getCytGraphNode() != 0) {
                        stageThreeCG(settings);
                } else if (settings.noCorrect()) {
                        stageThreeCompMult(settings);
                } else {
                        stageThreeCorrectGraph(settings);
                }
        } catch (exception& e) {
                cerr << "Fatal error: " << e.what() << endl;
                return EXIT_FAILURE;
        }

        cout << "Exiting... bye!" << endl;
        return EXIT_SUCCESS;
}
