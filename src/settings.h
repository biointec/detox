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

#ifndef SETTINGS_H
#define SETTINGS_H

#include <cstdlib>
#include <string>

#include "util.h"

// ============================================================================
// SETTINGS CLASS
// ============================================================================

class Settings {

private:
        // input file arguments
        std::string graphFilename;
        std::string readFilename;

        // options
        bool noCorrectFlag;     // do not correct graph, only compute multiplicities
        bool useQualFlag;       // use per-base quality scores (q-mers)
        double removeCoverage;  // in stage 1, remove all nodes with coverage <= to this threshold and all arcs with coverage 0

        int numThreads;         // number of threads

        int abundanceMin;       // min abundance threshold as passed to BCALM 2

        int crfDepth;           // CRF subgraph depth
        int crfMargin;          // CRF margin
        double crfFlow;         // CRF flow conservation strength
        double crfMaxFact;      // CRF maximum factor size

        double mmCov;           // initial coverage estimate
        double mmErrCov;        // initial error coverage estimate
        double mmODF;           // overdispersion factor
        int mmComponents;       // number of components in the mixture model
        int mmPloidy;           // ploidy [1 for haploid; 2 for diploid]

        int emMaxIter;          // EM maximum number of iterations
        double emConvEps;       // EM convergence tolerance (epsilon)
        double emTrainSize;     // EM number of nodes/arcs to train the model

        int cytGraphNode;       // Cytoscape graph center node

        int phredBase;          // ASCII value corresponding to Phred score Q=0

        /**
         * Print usage to stdout
         */
        void printUsage() const;

        /**
         * Print information about the program to stdout
         */
        void printProgramVersion() const;

public:
        /**
         * Constructor that parses program arguments
         * @param argc Argument count
         * @param argv Actual arguments
         */
        Settings(int argc, char** argv);

        /**
         * Get the graph filename
         * @return The graph filename
         */
        const std::string& getGraphFilename() const {
                return graphFilename;
        }

        /**
         * Get the graph filename produced in stage 1
         * @return The graph filename
         */
        std::string getStage1GraphFilename() const {
                return graphFilename + ".st1";
        }

        /**
         * Get the graph filename produced in stage 3
         * @return The graph filename
         */
        std::string getStage3GraphFilename() const {
                return graphFilename + ".st3";
        }

        /**
         * Get the node model filename produced in stage 2
         * @return The graph filename
         */
        std::string getNodeModelFilename() const {
                return "nodemodel.st2";
        }

        /**
         * Get the edge model filename produced in stage 2
         * @return The graph filename
         */
        std::string getEdgeModelFilename() const {
                return "edgemodel.st2";
        }

        /**
         * Get the true node multiplicity filename
         * @return The true node multiplicity filename
         */
        std::string getTrueNodeMultFilename() const {
                return "truemult.node";
        }

        /**
         * Get the true edge multiplicity filename
         * @return The true edge multiplicity filename
         */
        std::string getTrueEdgeMultFilename() const {
                return "truemult.edge";
        }

        /**
         * Get the estimated node multiplicity filename
         * @return The estimated node multiplicity filename
         */
        std::string getEstNodeMultFilename() const {
                return "estmult.node";
        }

        /**
         * Get the estimated edge multiplicity filename
         * @return The estimated edge multiplicity filename
         */
        std::string getEstEdgeMultFilename() const {
                return "estmult.edge";
        }

        /**
         * Get the node list filename
         * @return The node list filename
         */
        std::string getNodeListFilename() const {
                return "nodes.list";
        }

        /**
         * Get the edge list filename
         * @return The edge list filename
         */
        std::string getEdgeListFilename() const {
                return "edges.list";
        }

        /**
         * Get the read filename
         * @return The read filename
         */
        const std::string& getReadFilename() const {
                return readFilename;
        }

        /**
         * Check whether or not graph needs to be corrected
         * @return true when graph should not be corrected
         */
        bool noCorrect() const {
                return noCorrectFlag;
        }

        /**
         * Use per-base quality scores
         * @return true for q-mer counts, false for k-mer counts
         */
        bool useQual() const {
                return useQualFlag;
        }

        /**
         * Check if low coverage nodes and zero-coverage arcs should be removed in stage 1
         * @return the cutoff value for coverage removal in stage 1, -1 if no coverage removal required
         */
        double preCorrect() const {
                return removeCoverage;
        }

        /**
         * Get the number of threads
         * @return The number of threads
         */
        int getNumThreads() const {
                return numThreads;
        }

        /**
         * Get the CRF subgraph depth
         * @return The CRF subgraph depth
         */
        int getCRFDepth() const {
                return crfDepth;
        }

        /**
         * Get the CRF margin
         * @return The CRF margin
         */
        int getCRFMargin() const {
                return crfMargin;
        }

        /**
         * Get the CRF flow conservation strength
         * @return The CRF flow conservation strength
         */
        double getCRFFlowStrength() const {
                return crfFlow;
        }

        /**
         * Get the CRF maximum factor size
         * @return The CRF maximum factor size
         */
        double getCRFMaxFactSize() const {
                return crfMaxFact;
        }

        /**
         * Get the average coverage (to initialize EM)
         * @return The average coverage (to initialize EM)
         */
        double getInitCov() const {
                return mmCov;
        }

        /**
         * Get the average sequencing error coverage (to initialize EM)
         * @return The average sequencing error coverage (to initialize EM)
         */
        double getInitErrCov() const {
                return mmErrCov;
        }

        /**
         * Get the overdispersion factor (to initialize EM)
         * @return The overdispersion factor (to initialize EM)
         */
        double getInitODF() const {
                return mmODF;
        }

        /**
         * Get the number of components in the mixture model
         * @return The number of components in the mixture model
         */
        int getNumComp() const {
                return mmComponents;
        }

        /**
         * Get the average sequencing error coverage (to initialize EM)
         * @return The average sequencing error coverage (to initialize EM)
         */
        int getPloidy() const {
                return mmPloidy;
        }

        int getAbundanceMin() const {
                return abundanceMin;
        }

        /**
         * Get the maximum number of EM iterations
         * @return The maximum number of EM iterations
         */
        int getEMMaxIter() const {
                return emMaxIter;
        }

        /**
         * Get the convergence tolerance (epsilon) in EM
         * @return The convergence tolerance (epsilon) in EM
         */
        double getEMConvEps() const {
                return emConvEps;
        }

        /**
         * Get the number of nodes/arcs to train the mixture model using EM
         * @return The number of nodes/arcs to train the mixture model using EM
         */
        double getEMTrainingSize() const {
                return emTrainSize;
        }

        /**
         * Get the IO block size during multithreading
         * @return The IO block size during multithreading
         */
        size_t getThreadIOWorkSize() const {
                return 4096;
        }

        /**
         * Get the graph block size during multithreading
         * @return The graph block size during multithreading
         */
        size_t getThreadGraphWorkSize() const {
                return 1024;
        }

        /**
         * Get the Cytoscape graph center node (0 = none)
         * @return The Cytoscape graph center node
         */
        int getCytGraphNode() const {
                return cytGraphNode;
        }

        /**
         * Get the ASCII value corresponding to Phred score Q = 0
         * @return ASCII value corresponding to Phred score Q = 0
         */
        int getPhredBase() const {
                return phredBase;
        }

        /**
         * Check if necessary to execute stage 1
         * @return true or false
         */
        bool stageOneNecessary() const {
                return !Util::fileExists(getStage1GraphFilename());
        }

        /**
         * Check if necessary to execute stage 2
         * @return true or false
         */
        bool stageTwoNecessary() const {
                return !(Util::fileExists(getNodeModelFilename()) &&
                         Util::fileExists(getEdgeModelFilename()));
        }
};

#endif
