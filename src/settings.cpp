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

#include "settings.h"
#include "util.h"

#include <thread>
#include <iostream>
#include <fstream>

using namespace std;

void Settings::printProgramVersion() const
{
        cout << "Detox, a de Bruijn graph based assembler -- version "
             << DETOX_MAJOR_VERSION << "." << DETOX_MINOR_VERSION
             << "." << DETOX_PATCH_LEVEL << "\n";

        cout << "Copyright (C) 2019-2026 Jan Fostier and Aranka Steyaert\n"
                "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n"
                "This is free software; see the source for copying conditions. "
                "There is NO\nwarranty; not even for MERCHANTABILITY or "
                "FITNESS FOR A PARTICULAR PURPOSE." << endl;
}

void Settings::printUsage() const
{
        cout <<

        "Usage: detox [options] readfile\n\n  File "
        "  File \"readfile\" contains a list of input FASTQ files\n\n"

        " [options without arg]\n"
        "  -no-qual\t\tDo NOT use Phred quality scores to weigh k-mer coverage.\n\n"

        "  -help\t\t\tDisplay help page\n"
        "  -version\t\tDisplay version\n\n"

        " [options with 1 arg]\n"
        "  -k\t\t\tk-mer size [default = 21]\n"
        "  -num-threads\t\tNumber of threads [default = #cores]\n\n"

        "  -nn-model\t\tNeural network for multiplicity inference [default = \"\" (use CRF)]\n\n"

        "  -crf-nb-size\t\tCRF neighborhood size [default = 5]\n"
        "  -crf-margin\t\tCRF factor margin [default = 2]\n"
        "  -crf-flow\t\tCRF flow conservation strength [default = 1.0e7]\n"
        "  -crf-max-fact\t\tCRF maximum factor size [default = 1.0e6]\n\n"

        "  -mm-coverage\t\tInitial coverage est. in mixture model [default = auto]\n"
        "  -mm-err-cov\t\tInitial error coverage est. in mixture model [default = 1.0]\n"
        "  -mm-odf\t\tInitial overdispersion factor est. in mixture model [default = 1.5]\n"
        "  -mm-components\tNumber of components in mixture model [default = 6]\n\n"

        "  -em-max-iter\t\tMaximum number of EM iterations [default = 25]\n"
        "  -em-conv-eps\t\tRelative EM convergence epsilon [default = 1.0e-3]\n"
        "  -em-train-size\tNumber of nodes/arcs to use for EM training [default = 1e4]\n\n"

        "  -max-confl-degree\tMaximum conflict degree for repeat resolution [default = 0.25]\n"
        "  -min-count\t\tPaired-end read minimum support [default = 2]\n\n"

        "  -phred-base\t\tASCII value corresponding to Phred score Q = 0 [default = 33]\n\n"

        "Report bugs to Jan Fostier <jan.fostier@ugent.be>" << endl;
}

Settings::Settings(int argc, char ** argv) :
        useQualFlag(true), k(21), numThreads(thread::hardware_concurrency()),
        crfDepth(5), crfMargin(2), crfFlow(1e7), crfMaxFact(1e6), mmCov(-1.0),
        mmErrCov(1.0), mmODF(1.5), mmComponents(6), abundanceMin(1), emMaxIter(25),
        emConvEps(1e-3), emTrainSize(1e4), MCD(0.25), minCount(2), phredBase(33)
{
        const int reqArguments = 1;     // not counting argument 0

        // no arguments provided
        if (argc <= 1) {
                printUsage();
                exit(EXIT_FAILURE);
        }

        // only one argument provided
        if (argc == 1) {
                string arg(argv[1]);

                if (arg == "-help") {
                        printUsage();
                        exit(EXIT_SUCCESS);
                } else if (arg == "-version") {
                        printProgramVersion();
                        exit(EXIT_SUCCESS);
                }
        }

        // process optional arguments (we know there are at least two arguments)
        for (int i = 1; i < argc-reqArguments; i++) {
                string arg(argv[i]);

                // process options without arguments
                if (arg == "-no-qual") {
                        useQualFlag = false;
                        continue;
                } else if (arg == "-help") {
                        printUsage();
                        exit(EXIT_SUCCESS);
                } else if (arg == "-version") {
                        printProgramVersion();
                        exit(EXIT_SUCCESS);
                } else if (i == argc-reqArguments-1) {
                        cerr << "Unknown option or missing argument: "
                             << argv[i] << endl;
                        printUsage();
                        exit(EXIT_FAILURE);
                }

                // process options with arguments
                if (arg == "-k") {
                        numThreads = atoi(argv[i+1]);
                } else if (arg == "-num-threads") {
                        numThreads = atoi(argv[i+1]);
                } else if (arg == "-nn-model") {
                        nnModel = argv[i+1];
                } else if (arg == "-crf-nb-size") {
                        crfDepth = atoi(argv[i+1]);
                } else if (arg == "-crf-margin") {
                        crfMargin = atoi(argv[i+1]);
                } else if (arg == "-crf-flow") {
                        crfFlow = atof(argv[i+1]);
                } else if (arg == "-crf-max-fact") {
                        crfMaxFact = atof(argv[i+1]);
                } else if (arg == "-mm-coverage") {
                        mmCov = atof(argv[i+1]);
                } else if (arg == "-mm-err-cov") {
                        mmErrCov = atof(argv[i+1]);
                } else if (arg == "-mm-odf") {
                        mmODF = atof(argv[i+1]);
                } else if (arg == "-mm-components") {
                        mmComponents = atoi(argv[i+1]);
                } else if (arg == "-em-max-iter") {
                        emMaxIter = atoi(argv[i+1]);
                } else if (arg == "-em-conv-eps") {
                        emConvEps = atof(argv[i+1]);
                } else if (arg == "-em-train-size") {
                        emTrainSize = atof(argv[i+1]);
                } else if (arg == "-max-confl-degree") {
                        MCD = atof(argv[i+1]);
                } else if (arg == "-min-count") {
                        minCount = atoi(argv[i+1]);
                } else if (arg == "-phred-base") {
                        phredBase = atoi(argv[i+1]);
                } else {
                        cerr << "Unknown argument: " << argv[i] << endl;
                        printUsage();
                        exit(EXIT_FAILURE);
                }

                i++;    // if we reach this point, an argument was processed
        }

        readFilename = argv[argc-1];

        // check read file
        if (readFilename.empty()) {
                cerr << "Specify input read file" << endl;
                printUsage();
                exit(EXIT_FAILURE);
        }

        if (!Util::fileExists(readFilename))
                throw runtime_error("cannot open file " + readFilename);

        // perform sanity checks on input paramters
        if (k > 63)
                throw runtime_error("Maximum value for k is 63");
        if (k % 2 == 0)
                throw runtime_error("The k-mer size must be odd");
        if (numThreads < 1)
                throw runtime_error("Number of threads must be >= 1");
#ifndef HAVE_TORCH
        if (!nnModel.empty())
                throw runtime_error("Trying to pass a Torch model but Detox is built without Torch support. Please recompile.");
#endif
        if (crfDepth < 0 || crfDepth > 10)
                throw runtime_error("CRF subgraph depth must be "
                                    "in range [0..10]");
        if (crfMargin < 0 || crfMargin > 5)
                throw runtime_error("CRF factor margin must be "
                                    "in range [0..5]");
        if (crfFlow < 1.0)
                throw runtime_error("CRF flow conservation strength "
                                    "must be >= 1.0");
        if (crfMaxFact < 1e3)
                throw runtime_error("CRF maximum factor size must "
                                    "be >= 1e3");
        // Note: negative avgCov means "auto detect" and is hence allowed
        if (mmCov >= 0.0 && mmCov < 1.0)
                throw runtime_error("Mixture model initial coverage "
                                    "must be >= 1.0");
        if (mmCov >= 0.0 && mmCov <= mmErrCov)
                throw runtime_error("Mixture model initial coverage must be "
                                    "higher than the initial error coverage");
        if (mmErrCov <= 0.0)
                throw runtime_error("Mean sequencing error coverage "
                                    "must be > 0.0");
        if (mmODF < 1.0)
                throw runtime_error("Mixture model initial overdispersion "
                                    "factor must be >= 1.0");
        if (mmComponents < 2 || mmComponents > 15)
                throw runtime_error("Mixture model number of components must "
                                    "in range [2..15]");
        if (emMaxIter < 1)
                throw runtime_error("Maximum number of EM iterations "
                                    "must be >= 1");
        if (emConvEps >= 1.0 || emConvEps < 1e-8)
                throw runtime_error("EM relative converage tolerance "
                                    "should be in range [1e-8...1[");
        if (emTrainSize < 1e2)
                throw runtime_error("EM traing size should be >= 1e2");
        if (MCD < 0.0f || MCD > 1.0f)
                throw runtime_error("Maximum conflict degree should be in range [0..1]");
        if (minCount < 0)
                throw runtime_error("Paired-end read minimum support should be >= 0");
        if (phredBase < 0 || phredBase > 255)
                throw runtime_error("ASCII value for Phred score Q = 0 must "
                                    "be in range [0..255]");
}
