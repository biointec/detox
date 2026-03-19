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
#include "nninference.h"
#include "dbgraph.h"
#include "settings.h"

#include <iostream>

using namespace std;

// ============================================================================
// Node/Edge multiplicity estimation, either using CRF or using NN
// ============================================================================

std::tuple<NodeMap<int>, EdgeMap<int>>
MultEstimation::computeMultiplicitiesNN()
{
        NodeMap<int> estNodeMult;
        EdgeMap<int> estEdgeMult;

#ifdef HAVE_TORCH
        vector<NodeRep> nodeRep = dBG.getNodeReps(dBG.getNumValidNodes());
        vector<EdgeRep> edgeRep = dBG.getEdgeReps(dBG.getNumValidArcs());

        const float& nodeCov = nodeModelCRF.getLambda();
        const float& edgeCov = edgeModelCRF.getLambda();

        cout << "Infering multiplicities of " << nodeRep.size() << " nodes and "
             << edgeRep.size() << " edges" << endl;

        // create a node features tensor: reorder and expand features
        // such that they appear in node order
        // [-A, A, -B, B, -C, C, ... -N, N] with N = #nodes
        const int N = nodeRep.size();
        BiNodeMap<int> id2idx;
        torch::Tensor node_features = torch::empty({2*N, 3}, torch::kFloat32);
        for (size_t i = 0; i < N; i++) {
                NodeRep nr = nodeRep[i];
                SSNode n = dBG.getSSNode(nr);

                id2idx[-NodeID(nr)] = 2*i;
                id2idx[ NodeID(nr)] = 2*i+1;

                // compute average CG percentage
                string s = n.getSequence();
                int nCG = 0;
                for (const char c : s)
                        if (c == 'C' || c == 'G')
                                nCG++;
                float CGperc = (float)nCG / (float)s.size();

                node_features[2*i  ][0] = n.getAvgCov() / nodeCov; // normalize
                node_features[2*i+1][0] = n.getAvgCov() / nodeCov; // normalize
                node_features[2*i  ][1] = CGperc;
                node_features[2*i+1][1] = CGperc;
                node_features[2*i  ][2] = log(n.getMarginalLength());
                node_features[2*i+1][2] = log(n.getMarginalLength());
        }

        // create a {2, 2E} edge index tensor. For every edge (src, dst),
        // the reverse-complementary edge (-dst, -src) is added.
        const int E = edgeRep.size();
        torch::Tensor edge_index = torch::empty({2, 2*E}, torch::kInt32);
        for (int i = 0; i < E; i++) {
                NodeID src = edgeRep[i].getSrcID();
                NodeID dst = edgeRep[i].getDstID();

                edge_index[0][2*i  ] = id2idx[src]; // fwd edge
                edge_index[1][2*i  ] = id2idx[dst]; // fwd edge
                edge_index[0][2*i+1] = id2idx[-dst]; // rc edge
                edge_index[1][2*i+1] = id2idx[-src]; // rc edge
        }

        // create an edge features tensor
        torch::Tensor edge_features = torch::empty({2*E, 2}, torch::kFloat32);
        for (size_t i = 0; i < E; i++) {
                Arc& a = dBG.getArc(edgeRep[i]);
                SSNode src = dBG.getSSNode(edgeRep[i].getSrcID());
                string s = src.getRightKmer().str();
                SSNode dst = dBG.getSSNode(edgeRep[i].getDstID());
                char c = dst.peekNucleotideMarginalLeft();
                s += c;
                int nCG = 0;
                for (const char c : s)
                        if (c == 'C' || c == 'G')
                                nCG++;
                float CGperc = (float)nCG / (float)s.size();

                edge_features[2*i  ][0] = a.getCov() / edgeCov;
                edge_features[2*i+1][0] = a.getCov() / edgeCov;
                edge_features[2*i  ][1] = CGperc;
                edge_features[2*i+1][1] = CGperc;
        }

        // call forward on the model to do the actual inference
        vector<torch::jit::IValue> inputs;
        inputs.push_back(node_features);
        inputs.push_back(edge_index);
        inputs.push_back(edge_features);

        at::Tensor node_out, edge_out;
        try {
                auto output = modelNN.forward(inputs).toTuple();
                node_out = output->elements()[0].toTensor();
                edge_out = output->elements()[1].toTensor();
        } catch (const c10::Error& e) {
                throw std::runtime_error(std::string("Error during inference: ") + e.what());
        }

        // softmax + take class with highest probability
        node_out = torch::softmax(node_out, 1).argmax(1);
        edge_out = torch::softmax(edge_out, 1).argmax(1);

        for (size_t i = 0; i < N; i++)
                estNodeMult[nodeRep[i]] = node_out[2*i].item<int>();
        for (size_t i = 0; i < E; i++)
                estEdgeMult[edgeRep[i]] = edge_out[2*i].item<int>();

        // recompute multiplicity == maxMult using CRFs
        NodeMap<Multiplicity> nodeMult_m;
        for (auto [nr, m] : estNodeMult)
                if (m >= maxMult)
                        nodeMult_m[nr] = Multiplicity();

        EdgeMap<Multiplicity> edgeMult_m;
        for (auto [er, m] : estEdgeMult)
                if (m >= maxMult)
                        edgeMult_m[er] = Multiplicity();

        myCRFMult.computeMult(nodeMult_m, nodeModelCRF,
                              edgeMult_m, edgeModelCRF);

        for (const auto& [nr, mp] : nodeMult_m)
                estNodeMult[nr] = mp.getExpMult();
        for (const auto& [er, mp] : edgeMult_m)
                estEdgeMult[er] = mp.getExpMult();

#endif

        return {estNodeMult, estEdgeMult};
}

std::tuple<NodeMap<int>, EdgeMap<int>>
MultEstimation::computeMultiplicitiesCRF()
{
        vector<NodeRep> nodeRep = dBG.getNodeReps(dBG.getNumValidNodes());
        vector<EdgeRep> edgeRep = dBG.getEdgeReps(dBG.getNumValidArcs());

        // create multiplicity maps
        NodeMap<Multiplicity> nodeMult_m;
        for (size_t i = 0; i < nodeRep.size(); i++)
                nodeMult_m[nodeRep[i]] = Multiplicity();

        EdgeMap<Multiplicity> edgeMult_m;
        for (size_t i = 0; i < edgeRep.size(); i++)
                edgeMult_m[edgeRep[i]] = Multiplicity();

        myCRFMult.computeMult(nodeMult_m, nodeModelCRF,
                              edgeMult_m, edgeModelCRF);

        // convert multiplicities to integers
        NodeMap<int> nodeMult = NodeMap<int>(nodeMult_m.size());
        for (const auto& [id, m] : nodeMult_m)
                nodeMult[id] = m.getExpMult();
        EdgeMap<int> edgeMult = EdgeMap<int>(edgeMult_m.size());
        for (const auto& [id, m] : edgeMult_m)
                edgeMult[id] = m.getExpMult();

        return {nodeMult, edgeMult};
}

MultEstimation::MultEstimation(const DBGraph& dBG, const Settings& settings,
                               const string& nodeModelFilenameCRF,
                               const string& edgeModelFilenameCRF,
                               const string& modelFilenameNN) :
        dBG(dBG), haveNN(false),
        nodeModelCRF(nodeModelFilenameCRF), edgeModelCRF(edgeModelFilenameCRF),
        myCRFMult(dBG, settings.getCRFDepth(), settings.getCRFMargin(),
                  settings.getCRFMaxFactSize(), settings.getCRFFlowStrength(),
                  settings.getNumThreads(), settings.getCRFWorkSize())
{
        // try and load the NN model, if provided
        if (modelFilenameNN.empty())
                return;

#ifdef HAVE_TORCH
        try {
                modelNN = torch::jit::load(modelFilenameNN);
                modelNN.eval();
                cout << "Loaded Torch model: " << modelFilenameNN << endl;
                haveNN = true;
        } catch (const c10::Error& e) {
                throw runtime_error(string("Error loading the model: ") + e.what());
        }
#endif
}

std::tuple<NodeMap<int>, EdgeMap<int>>
MultEstimation::computeMultiplicities()
{
        return (haveNN) ?
                computeMultiplicitiesNN() :
                computeMultiplicitiesCRF();
}
