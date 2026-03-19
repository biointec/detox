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
#ifndef NNINFERENCE_H
#define NNINFERENCE_H

#ifdef HAVE_TORCH
#include <torch/script.h> // TorchScript header
#endif
#include <tuple>

#include "coverage.h"
#include "ssnode.h"
#include "crfmult.h"

class DBGraph;
class Settings;

// ============================================================================
// Node/Edge multiplicity estimation, either using CRF or using NN
// ============================================================================

class MultEstimation
{
private:
        const int maxMult = 5;
        const DBGraph& dBG;     // de Bruijn graph
        bool haveNN;            // true if succesfully loaded a NN model
#ifdef HAVE_TORCH
        torch::jit::script::Module modelNN;
#endif  // HAVE_TORCH
        CovModel nodeModelCRF;  // CRF node model
        CovModel edgeModelCRF;  // CRF edge model
        CRFMult myCRFMult;      // CRF multiplicity inference

        /**
         * @brief Compute all node/edge multiplicities for a de Bruijn graph
         * using the neural network model
         * @return tuple [NodeMap, EdgeMap] containing the multiplicities
         */
        std::tuple<NodeMap<int>, EdgeMap<int>>
                computeMultiplicitiesNN();

        /**
         * @brief Compute all node/edge multiplicities for a de Bruijn graph
         * using the CRF model
         * @return tuple [NodeMap, EdgeMap] containing the multiplicities
         */
        std::tuple<NodeMap<int>, EdgeMap<int>>
                computeMultiplicitiesCRF();

public:
        /**
         * @brief Constructor
         * @param dBG Const-reference to de Bruijn graph
         * @param settings Settings object
         * @param nodeModelFilenameCRF Filename to CRF node model
         * @param edgeModelFilenameCRF Filename to CRF edge model
         * @param modelFilenameNN Filename to NN model (empty if not used)
         */
        MultEstimation(const DBGraph& dBG, const Settings& settings,
                       const std::string& nodeModelFilenameCRF,
                       const std::string& edgeModelFilenameCRF,
                       const std::string& modelFilenameNN);

        /**
         * Compute all node/edge multiplicities for a de Bruijn graph
         * @return tuple [NodeMap, EdgeMap] containing the multiplicities
         */
        std::tuple<NodeMap<int>, EdgeMap<int>>
                computeMultiplicities();
};

#endif // NNINFERENCE_H
