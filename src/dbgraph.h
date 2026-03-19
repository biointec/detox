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

#ifndef DBGRAPH_H
#define DBGRAPH_H

#include "kmer/tkmer.h"
#include "ssnode.h"
#include "nodechain.h"
#include "kmernpp.h"
#include "library.h"
#include "readfile/fastq.h"
#include "readfile/fasta.h"

#include <mutex>
#include <vector>
#include <queue>

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class Settings;
class LibraryContainer;
class Multiplicity;
class PathInfo;
class ChainColl;
class NodeChainSet;
class Contig;
class NodeTracker;
class CRFMult;
class CovModel;
class PathFinder;

// ============================================================================
// DIJKSTRA AUXILIARY CLASSES
// ============================================================================

typedef std::pair<NodeID, int> NodeDepth;

struct NodeDepthComp {

        /**
         * Compare by depth (use greater because priority_queue.top() returns
         * by default the greatest element. We want to return the smallest.)
         */
        bool operator()(const NodeDepth& f, const NodeDepth& s) {
                return f.second > s.second;
        }
};

// ============================================================================
// PARALLEL PATH AUXILIARY CLASS
// ============================================================================

class ParallelPath {

public:
        std::vector<NodeID> orig;       // original path
        std::vector<NodeID> para;       // parallel path
        size_t alnOrig;                 // lenght of original path alignment
        size_t alnPara;                 // length of parallel path alignment

        ParallelPath() : alnOrig(0), alnPara(0) {};

        ParallelPath(std::vector<NodeID> orig, std::vector<NodeID> para,
                     size_t alnOrig, size_t alnPara) : orig(orig), para(para),
                     alnOrig(alnOrig), alnPara(alnPara) {}
};

// ============================================================================
// DE BRUIJN GRAPH CLASS
// ============================================================================

class DBGraph {

private:
        // ====================================================================
        // VARIABLES
        // ====================================================================

        const Settings& settings;       // reference to settings object

        DSNode *nodes;                  // graph nodes
        Arc *arcs;                      // graph arcs

        NodeID numNodes;                // number of nodes
        ArcID numArcs;                  // number of arcs

        NodeID numValidNodes;           // number of nodes
        ArcID numValidArcs;             // number of arcs

        // ====================================================================
        // GRAPH MODIFICATION
        // ====================================================================

        /**
         * Remove a node from the graph
         * @param nodeID node identifier
         * @param concatentate Concatenate adjacent nodes
         **/
        void removeNode(NodeID nodeID, bool concatenate = false);

public:
        void createNodeMapping(const DBGraph& sg, NodeID srcID, NodeID newSrcID,
                               BiNodeMap<NodePosPair>& mapping) const;

        /**
         * Remove an arc between two nodes
         * @param leftID Identifier of the left node
         * @param rightID Identifier of the right node
         * @param concatenate Concatenate adjacent nodes
         */
        void removeArc(NodeID leftID, NodeID rightID, bool concatenate = false);

        /**
         * Remove all right arcs from a node
         * @param nodeID Identifier of the node
         */
        void removeRightArcs(NodeID nodeID);

        /**
         * Remove all left arcs from a node
         * @param nodeID Identifier of the node
         */
        void removeLeftArcs(NodeID nodeID);

        /**
         * Concatenation the linear path around a seed node
         * @param seedID Node identifier of the seed node
         * @param nodeListv Linear path that was concatenated
         */
        void concatenateAroundNode(NodeID seedID,
                                   std::vector<NodeID>& nodeListv);

        // ====================================================================
        // COVERAGE
        // ====================================================================

        /**
         * Count node and arc coverage for a single read record
         * @param rr Read record
         * @param table < Kmer, NPP > table
         */
        void covCount(const FastQRecord& rr,
                      const KmerNPPTable& table);

        /**
         * Thread entry point for counting node and arc coverage
         * @param inputs Library container with input read files
         * @param table < Kmer, NPP > table with all contig endpoint kmers
         */
        void covCountThread(FastQReader& inputs,
                            const KmerNPPTable& table);

public:
        /**
         * Default constructor
         * @param settings Const-ref to settings object
         */
        DBGraph(const Settings& settings) : settings(settings), nodes(NULL),
                arcs(NULL), numNodes(0), numArcs(0) {};

        /**
         * Given a dBG and a subset of nodes and edges, build a subgraph
         * @param dBG Parent de Bruijn graph
         * @param subNodes Subset of nodes
         * @param subEdges Subset of edges
         * @param subToOrig Mapping nodes new to original
         * @param origToSub Mapping nodes original to new
         */
        DBGraph(const DBGraph& dBG,
                const NodeMap<int>& subNodes, const EdgeMap<int>& subEdges,
                BiNodeMap<NodeID>& subToOrig, BiNodeMap<NodeID>& origToSub);

        /**
         * Destructor
         */
        ~DBGraph() {
                delete [] nodes;
                delete [] arcs;
        }

        // ====================================================================
        // UTILITY ROUTINES
        // ====================================================================

        /**
         * Get the abundanceMin value from settings object
         * @return First coverage value in the histogram to use for fitting
         */
        int getAbundanceMin() const ;

        /**
         * Get the number of nodes
         * @return The number of nodes
         */
        NodeID getNumNodes() const {
                return numNodes;
        }

        /**
         * Get the number of arcs
         * @return The number of arcs
         */
        ArcID getNumArcs() const {
                return numArcs;
        }

        /**
         * Get the number of valid nodes
         * @return The number of valid nodes
         */
        NodeID getNumValidNodes() const {
                return numValidNodes;
        }

        /**
         * Get the number of valid arcs
         * @return The number of valid arcs
         */
        ArcID getNumValidArcs() const {
                return numValidArcs;
        }

        /**
         * Get the arc identifier given the arc pointer
         * @param arc Arc pointer
         * @return The arc identifier
         */
        size_t getArcID(Arc *arc) const {
                assert((arc - arcs) > 0);
                assert((arc - arcs) <= (std::ptrdiff_t)numArcs);

                return arc - arcs;
        }

        /**
         * Get the arc identifier given the edge representative
         * @param e Edge representative
         * @return The arc identifier
         */
        ArcID getArcID(const EdgeRep& e) const {
                Arc* arc = getSSNode(e.getSrcID()).rightArc(e.getDstID());
                return getArcID(arc);
        }

        /**
         * Get a reference to an arc, given the edge representative
         * @param e Edge representative
         * @return Reference to the arc
         */
        Arc& getArc(const EdgeRep& e) const {
                return *getSSNode(e.getSrcID()).rightArc(e.getDstID());
        }

        /**
         * Get a reference to an arc, given source and destination node
         * @param srcID Source node identifier
         * @param dstID Destination node identifier
         * @return Reference to the arc
         */
        Arc& getArc(NodeID srcID, NodeID dstID) const {
                return *getSSNode(srcID).rightArc(dstID);
        }

        /**
         * Get the coverage of an arc
         * @param srcID Source node identifier
         * @param dstID Destination node identifier
         * @return Arc coverage
         */
        Coverage getArcCov(NodeID srcID, NodeID dstID) const {
                return getArc(srcID,  dstID).getCov();
        }

        /**
         * Set the coverage of an arc
         * @param srcID Source node identifier
         * @param dstID Destination node identifier
         * @param cov Target coverage
         */
        void setArcCov(NodeID srcID, NodeID dstID, Coverage cov) {
                getArc( srcID,  dstID).setCov(cov);
                getArc(-dstID, -srcID).setCov(cov);
        }

        /**
         * Increase the coverage of an arc
         * @param srcID Source node identifier
         * @param dstID Destination node identifier
         * @param addCov Coverage to add
         */
        void incArcCov(NodeID srcID, NodeID dstID, Coverage addCov) {
                setArcCov(srcID, dstID,
                          getArc(srcID, dstID).getCov() + addCov);
        }

        /**
         * Check whether a node exists
         * @param nr Node representative
         * @return true or false
         */
        bool nodeExists(NodeRep nr) const {
                if (nr.getNodeID() == 0)
                        return false;
                if (abs(nr.getNodeID()) > getNumNodes())
                        return false;
                return getSSNode(nr.getNodeID()).isValid();
        }

        /**
         * Check whether a edge exists (and its adjacent nodes)
         * @param er Edge representative
         * @return true of false
         */
        bool edgeExists(EdgeRep er) const {
                if (!nodeExists(er.getSrcID()))
                        return false;
                if (!nodeExists(er.getDstID()))
                        return false;
                return getSSNode(er.getSrcID()).rightArc(er.getDstID()) != NULL;
        }

        /**
         * Check whether a path exists in the graph
         * @param nc Node chain
         * @return true or false
         */
        bool pathExists(const NodeChain& nc) const {
                for (size_t i = 0; i < nc.size(); i++)
                        if (!nodeExists(nc[i]))
                                return false;
                for (size_t i = 1; i < nc.size(); i++)
                        if (getSSNode(nc[i-1]).rightArc(nc[i]) == NULL)
                                return false;
                return true;
        }

        /**
         * Set the flag1 of all nodes
         * @param value Target value
         */
        void setAllFlags1(bool value) {
                for (NodeID i = 1; i <= numNodes; i++)
                        getSSNode(i).setFlag1(value);
        }

        /**
         * Set the flag2 of all nodes
         * @param value Target value
         */
        void setAllFlags2(bool value) {
                for (NodeID i = 1; i <= numNodes; i++)
                        getSSNode(i).setFlag2(value);
        }

        /**
         * Given an arc identifier, that the adjacent nodes
         * @param arcID Arc identifier
         * @return pair of adjacent nodes < left, right >
         */
        std::pair<NodeID, NodeID> getAdjacentNodes(size_t arcID) const {
                const Arc& arc = getArc(arcID);
                NodeID rightID = arc.getNodeID();

                SSNode right = getSSNode(rightID);
                for (ArcIt it = right.leftBegin(); it != right.leftEnd(); it++) {
                        NodeID leftID = it->getNodeID();
                        SSNode left = getSSNode(leftID);

                        if (left.rightArc(rightID) == &arc)
                                return std::make_pair(leftID, rightID);
                }

                NodeID leftID = arc.getNodeID();
                SSNode left = getSSNode(leftID);
                for (ArcIt it = left.rightBegin(); it != left.rightEnd(); it++) {
                        rightID = it->getNodeID();
                        right = getSSNode(rightID);

                        if (right.leftArc(leftID) == &arc)
                                return std::make_pair(leftID, rightID);
                }

                return std::make_pair(0, 0);
        }

        /**
         * Given an arc identifier, get the edge representative
         * @param arcID Arc identifier
         * @return Edge representative
         */
        EdgeRep getEdgeRep(ArcID arcID) const {
                std::pair<NodeID, NodeID> nodes = getAdjacentNodes(arcID);
                return EdgeRep(nodes.first, nodes.second);
        }

        /**
         * Get a reference to an arc, given the arcID
         * @param arcID Identifier for the arc
         * @return Reference to the arc
         */
        Arc& getArc(size_t arcID) const {
                return arcs[arcID];
        }

        /**
         * Get a reference to a double stranded node, given the nodeID
         * @param nodeID Identifier for the node
         * @return Reference to the node
         */
        DSNode& getDSNode(NodeID nodeID) const {
                assert(nodeID > 0 && nodeID <= numNodes);
                return nodes[nodeID];
        }

        /**
         * Get a single stranded node, given the nodeID
         * @param nodeID Identifier for the node
         * @return A Single Stranded node
         */
        SSNode getSSNode(NodeID nodeID) const {
                NodeID uNodeID = abs(nodeID);
                assert(uNodeID != 0 && uNodeID <= numNodes);
                return SSNode(nodes + uNodeID, nodeID);
        }

        std::vector<std::pair<NodeID, double>> getRightNeigbors(NodeID id) const;

        /**
         * Check whether two NodePosPairs are consecutive in the graph
         * @param left Left NodePosPair
         * @param right Right NodePosPair
         * @param offset Distance between the pairs
         * @return True if the NodePosPairs are consecutive
         */
        bool consecutiveNPP(const NodePosPair& left,
                            const NodePosPair& right, size_t offset = 1) const;

        /**
         * Check whether two NodePosPairs cross an arc in the graph
         * @param left Left NodePosPair
         * @param right Right NodePosPair
         * @return True if the NodePosPairs cross an arc
         */
        bool crossesArc(const NodePosPair& left,
                        const NodePosPair& right) const
        {
                if (!consecutiveNPP(left, right))
                        return false;

                return ( (left.getPosition()+1) != right.getPosition() );
        }

        /**
         * Defragment the node array
         */
        void defragNodes();

        /**
         * Defragment the arc array
         */
        void defragArcs();

        /**
         * Populate a <kmer, NodePosPair> table
         * @param table Table to populate
         */
        void createKmerNPPTable(KmerNPPTable& table) const;

        // ====================================================================
        // I/O ROUTINES
        // ====================================================================

        /**
         * Load a graph from the BCALM 2 format
         * @param filename Input filename
         */
        void loadBCalm(const std::string& filename);

        void loadGraph(const std::string& filename);

        /**
         * Write a graph in the BCALM 2 format
         * @param filename Input filename
         */
        void writeBCalm(const std::string& filename);

        /**
         * Write graph in binary format
         * @param filename Output filename
         */
        void writeBinary(const std::string& filename) const;

        /**
         * Load graph in binary format
         * @param filename Input filename
         */
        void loadBinary(const std::string& filename);

        /**
         * Write the graph's contigs in FASTA format
         * @param filename Output filename
         * @param contigs Contigs
         */
        void writeContigs(const std::string& filename,
                          const std::vector<NodeChain>& contigs,
                          const std::vector<std::string>& contigNames) const;

        /**
         * Write a Cytoscape graph of the current graph
         * @param filename Filename of the cytoscape graph
         * @param nodes Nodes in subgraph
         * @param edges Edges in subgraph
         * @param estNodeMult Estimated node multiplicity
         * @param estEdgeMult Estimated edge multiplicity
         * @param trueNodeMult True node multiplicity
         * @param trueEdgeMult True edge multiplicity
         */
        void writeCytoscapeGraph(const std::string& filename,
                                 std::vector<NodeID> nodes,
                                 std::vector<EdgeID> edges,
                                 const NodeMap<Multiplicity>& estNodeMult,
                                 const EdgeMap<Multiplicity>& estEdgeMult,
                                 const NodeMap<int>& trueNodeMult,
                                 const EdgeMap<int>& trueEdgeMult) const;

        void writeCytoscapeSubgraph(const std::string& filename,
                                    std::vector<NodeID> nodes,
                                    std::vector<EdgeID> edges,
                                    BiNodeMap<NodeID> subToOrig) const;

        /**
         * Check the validity of the de Bruijn graph
         */
        void sanityCheck();

        /**
         * Get the entire de Bruijn graph. Output will contain both nodeID
         * and -nodeID.
         * @param nodes Vector of nodes in the graph (output)
         * @param edges Vector of edges in the graph (output)
         */
        void getGraph(std::vector<NodeID>& nodes, std::vector<EdgeID>& edges);

        /**
         * Get a directed subgraph of the de Bruijn graph. Output can contain
         * both nodeID and -nodeID in presence of palindromic repeats
         * @param seedNode Seed node
         * @param nodes Vector of nodes in the subgraph (output)
         * @param edges Vector of edges in the subgraph (output)
         * @param maxDepth Maximum depth (in terms of number of nodes)
         */
        void getSubgraphAux(
                std::priority_queue<NodeDepth, std::vector<NodeDepth>, NodeDepthComp>& todo,
                std::vector<NodeID>& nodes, std::vector<EdgeID>& edges,
                size_t maxDepth = 0) const;

        /**
         * Get a directed subgraph of the de Bruijn graph. Output can contain
         * both nodeID and -nodeID in presence of palindromic repeats
         * @param seedNode Seed node
         * @param nodes Vector of nodes in the subgraph (output)
         * @param edges Vector of edges in the subgraph (output)
         * @param maxDepth Maximum depth (in terms of number of nodes)
         */
        void getSubgraph(NodeID seedNode, std::vector<NodeID>& nodes,
                         std::vector<EdgeID>& edges, size_t maxDepth = 0) const;

        /**
         * Get a directed subgraph of the de Bruijn graph. Output can contain
         * both nodeID and -nodeID in presence of palindromic repeats.
         * Used for Cytoscape visualization
         * @param seedChain Seed chain
         * @param nodes Vector of nodes in the subgraph (output)
         * @param edges Vector of edges in the subgraph (output)
         * @param maxDepth Maximum depth (in terms of number of nodes)
         */
        void getSubgraph(const NodeChain& seedChain, std::vector<NodeID>& nodes,
                         std::vector<EdgeID>& edges, size_t maxDepth = 0) const;

        /**
         * Create a vector of node representations stored in a file
         * @param filename File name
         * @return A vector of node representations
         */
        std::vector<NodeRep> getNodeReps(const std::string& filename) const;

        /**
         * Create a vector of N randomly selected node representations
         * @param N The number of desired elements
         * @return A vector of node representations
         */
        std::vector<NodeRep> getNodeReps(size_t N) const;

        /**
         * Create a vector of edge representations stored in a file
         * @param filename File name
         * @return A vector of edge representations
         */
        std::vector<EdgeRep> getEdgeReps(const std::string& filename) const;

        /**
         * Create a vector of N randomly selected edge representations
         * @param N The number of desired elements
         * @return A vector of edge representations
         */
        std::vector<EdgeRep> getEdgeReps(size_t N) const;

        /**
         * Get a vector of node identifiers which are right-branching
         * @return A vector of node identifiers which are right-branching
         */
        std::vector<NodeID> getBranchingNodes() const {
                std::vector<NodeID> result;
                for (NodeID id = -numNodes; id <= numNodes; id++) {
                        if (!nodeExists(id))
                                continue;
                        if (getSSNode(id).numRightArcs() > 1)
                                result.push_back(id);
                }
                return result;
        }


        std::vector<NodeRep> getLowCovForks(double threshold) const;
        std::vector<NodeRep> getLowCovForks(double loCov, double hiCov) const;
        std::vector<NodeRep> getNonInterfering(double threshold);

        /**
         * Get nodes with an average coverage below the threshold that are
         * adjacent to high-coverage nodes
         * @param loCov Low coverage threshold
         * @param hiCov High coverage threshold
         * @return A vector of node representations
         */
        std::vector<NodeRep> getLowCovNodesAdjToHighCov(double loCov,
                                                        double hiCov) const;

        /**
         * Get edges with a coverage below the threshold
         * @param threshold Coverage threshold
         * @return A vector of edge representations
         */
        std::vector<EdgeRep> getLowCovEdges(double threshold) const;

        /**
         * Remove tips from graph. Coverage is transfered to a parallel path.
         * @param covCutoff Maximum coverage for the tip
         * @param maxLen Maximum marginal length
         * @param toKeep Set of nodes that should not be removed
         * @param extendTips Merge longer tips into shorter ones
         * @return Number of nodes removed
         */
        size_t removeLowCovTips(double covCutoff, size_t maxLen,
                                const std::set<NodeRep>& toKeep = {},
                                bool extendTips = true);

        /**
         * Remove bubbles from graph and concatenate
         * @param covCutoff Maximum coverage for the bubble
         * @param maxLen Maximum marginal length
         * @param toKeep Set of nodes that should not be removed
         * @return Number of nodes removed
         */
        size_t removeLowCovBubbles(double covCutoff, size_t maxLen,
                                   const std::set<NodeRep>& toKeep = {});

        bool getTipsTree(NodeID srcID, std::set<NodeID>& tips, size_t maxLen);

        /**
         * Remove tips from graph.
         * @param maxLen Maximum marginal length
         * @return Number of tips removed
         */
        size_t removeTipsCRF(size_t maxLen, CRFMult& myCRFMult,
                             CovModel& nodeModel, CovModel& edgeModel);

        /**
         * Remove tips from graph.
         * @param maxLen Maximum marginal length
         * @return Number of tips removed
         */
        size_t removeBubblesCRF(size_t maxLen, CRFMult& myCRFMult,
                                CovModel& nodeModel, CovModel& edgeModel);

        /**
         * Remove tips from graph.
         * @param maxLen Maximum marginal length
         * @return Number of tips removed
         */
        size_t removeNodesCRF(size_t maxLen, CRFMult& myCRFMult,
                              CovModel& nodeModel, CovModel& edgeModel);

        /**
         * Get nodes with an average coverage below the threshold
         * @param threshold Coverage threshold
         * @return A vector of node representations
         */
        std::vector<NodeRep> getLowCovNodes(double threshold) const;

        void filterNodes(const std::vector<NodeRep>& nodeReps,
                         const std::vector<bool>& flowOK,
                         std::vector<NodeRep>& toRemove,
                         std::vector<NodeRep>& toKeep,
                         int depth);

        /**
         * Filter non-interfering nodes
         * @param nodeMult The multiplicity of a set of nodes (input)
         * @param toRemove The nodes that can be removed
         * @param toKeep The nodes in the neighborhood of those that are removed
         * @param depth Subgraph depth
         */
        void filterNodes(const NodeMap<Multiplicity>& nodeMult,
                         std::vector<NodeRep>& toRemove,
                         std::vector<NodeRep>& toKeep,
                         int depth);

        void filterEdges(const std::vector<EdgeRep>& edgeReps,
                         const std::vector<bool>& flowOK,
                         std::vector<EdgeRep>& toRemove,
                         std::vector<EdgeRep>& toKeep,
                         int depth);

        /**
         * Remove nodes with expected coverage == 0 and conservation-of-flow
         * @param covCutoff Maximum coverage for the node
         */
        void removeLowCovNodesFlow(double covCutoff, CRFMult& myCRFMult,
                                   CovModel& nodeModel, CovModel& edgeModel);

        /**
         * Remove edges with expected coverage == 0 and conservation-of-flow
         * @param covCutoff Maximum coverage for the node
         */
        void removeLowCovEdgesFlow(double covCutoff, CRFMult& myCRFMult,
                                   CovModel& nodeModel, CovModel& edgeModel);

        /**
         * Remove nodes with CRF coverage == 0
         * @param covCutoff Maximum coverage for the node
         */
        void removeLowCovNodesCRF(double covCutoff, CRFMult& myCRFMult,
                                  CovModel& nodeModel, CovModel& edgeModel);


        /**
         * Filter non-interfering edges
         * @param edgeMult The multiplicity of a set of edges (input)
         * @param toRemove The edges that can be removed
         * @param toKeep The edges in the neighborhood of those that are removed
         * @param depth Subgraph depth
         */
        void filterEdges(const EdgeMap<Multiplicity>& edgeMult,
                         std::vector<EdgeRep>& toRemove,
                         std::vector<EdgeRep>& toKeep,
                         int depth);

        /**
         * Remove edges with CRF coverage == 0
         * @param covCutoff Maximum coverage for the edge
         */
        void removeLowCovEdgesCRF(double covCutoff, CRFMult& myCRFMult,
                                  CovModel& nodeModel, CovModel& edgeModel);

        /**
         * Get chimeric connections with an average coverage below the threshold
         * @param threshold Coverage threshold
         * @param maxLen Maximum marginal length
         * @return A vector of node representations
         */
        std::vector<NodeRep> getLowCovChimConn(double threshold,
                                               size_t maxLen) const;

        void smoothBubbles(size_t maxLen, size_t maxED);

        double getInitialKmerCovEstimate(double errLambda, double p) const;

        /**
         * Calculate the coverage
         * @param libraries Input libraries (= read files)
         * @param table Populated <kmer, nodePosPair> table
         */
        void getCovFromReads(LibraryContainer& libraries, const KmerNPPTable& table);

        /**
         * Invalidate low-coverage nodes and arcs from a graph
         * @param covCutoff Maximum coverage of a node/arc to delete
         * @param maxMargLength Maximum marginal length of a node to delete
         */
        void removeCoverage(double covCutoff, size_t maxMargLength);

        /**
         * Invalidate low-coverage nodes (and zero coverage arcs) from a graph
         * @param covCutoff Maximum coverage of a node to delete
         * @param maxMargLength Maximum marginal length of a node to delete
         */
        void removeCoverageNodes(double covCutoff, size_t maxMargLength);

        /**
         * Convert a node chain to a string
         * @param nc A node chain
         * @return An stl string (output)
         */
        std::string convertNodesToString(const NodeChain& nc) const;

        /**
         * Convert a node chain to a string: only the last k-mer of the first
         * node and the first k-mer of the last node are included
         * @param nc A node chain
         * @return An stl string (output)
         */
        std::string convertNodesToShortestString(const NodeChain& nc) const;

        /**
         * Convert a shortest string to a node chain
         * @param srcID Node identifier that ends with first k-mer of str
         * @param str An stl string
         * @return A node chain (output)
         */
        NodeChain convertShortestStringToNodes(NodeID srcID,
                const std::string& str) const;

        ConsNodeChain convertCNC(const ConsNodeChain& srcCNC,
                                 const DBGraph& newGraph,
                                 const BiNodeMap<NodeID>& thisToNew) const;

        NodeChain convertNC(const NodeChain& thisNC,
                            const DBGraph& newGraph,
                            const BiNodeMap<NodePosPair>& origToSub) const;

        /**
         * Concatentate linear paths
         * @return True if at least one node was merged
         **/
        bool concatenateNodes(const std::set<NodeRep>& toKeep = {});

        /**
         * Remove of list of nodes from the graph
         * @param nodes Vector of node representatives to remove
         * @param concatenate Concatenate adjacent nodes?
         */
        void removeNodes(const std::vector<NodeRep>& nodes,
                         bool concatenate = false);

        /**
         * Try to glue tips together
         * @param libraries Input libraries (= read files)
         */
        void glueTips(LibraryContainer& libraries);

        /**
         * Glue two tips
         * @param leftID Identifier of the left tip
         * @param rightID Identifier of the right top
         */
        void glueTips(NodeID leftID, NodeID rightID);

        /**
         * Remove of list of edges from the graph
         * @param nodes Vector of edge representatives to remove
         * @param concatenate Concatenate adjacent nodes?
         */
        void removeEdges(const std::vector<EdgeRep>& edges,
                         bool concatenate = false);

        /**
         * @brief Update chain collection by adding a single (isolated) read
         * @param read Read to add
         * @param noi Nodes-of-interest (= nodes for which read info is kept)
         * @param chainColl Read collection (updated)
         */
        void addReadInfo(const NodeChain& read,
                         const BiNodeMap<bool>& noi,
                         BiNodeMap<ChainColl>& chainColl) const;

        /**
         * @brief Update read information by adding a library container
         * @param library Read library
         * @param noi Nodes-of-interest (= nodes for which read info is kept)
         * @param readInfo Read information data structure (updated)
         */
        void addReadInfo(const LibraryContainer& library,
                         const BiNodeMap<bool>& noi,
                         BiNodeMap<PathInfo>& readInfo) const;

        /**
         * @brief Update read information by adding a paired-end read
         * @param L First read
         * @param R Second read
         * @param readID Read identifier
         * @param noi Nodes-of-interest (= nodes for which read info is kept)
         * @param readInfo Read information data structure (updated)
         */
        void addPERInfo(const NodeChain& L, const NodeChain& R,
                        size_t readID, const BiNodeMap<bool>& noi,
                        BiNodeMap<PathInfo>& readInfo) const;

        /**
         * @brief Update read information by adding a library container
         * @param noi Nodes-of-interest (= nodes for which read info is kept)
         * @param readInfo Read information data structure (updated)
         */
        void addPERInfo(const LibraryContainer& library,
                        const BiNodeMap<bool>& noi,
                        BiNodeMap<PathInfo>& readInfo) const;

        /**
         * @brief Build paired-end read information for a collection of contigs
         * @param contigs Contigs of interest
         * @param node2contig Map that links (unique) nodes to their contig
         * @param readInfo Read information
         * @param MCD Maximum conflict degree
         * @param perInfo Paired-end read information (output)
         */
        void buildPerInfo(const LibraryContainer& library,
                          const BiNodeMap<Contig>& contigs,
                          const BiNodeMap<int>& node2contig,
                          BiNodeMap<PathInfo>& readInfo, float MCD,
                          NodeChainSet& u2uPer) const;

        /**
         * @brief Filter path information: 1) filter direct paths for which the
         * conflict degree <= MCD and 2) filter paired-end reads <= minCount
         * @param pathInfo Path information (input / output)
         * @param MCD Maximum conflict degree
         * @param minCount Minimum count
         */
        void filterPathInfo(BiNodeMap<PathInfo>& pathInfo,
                            float MCD, int minCount) const;

        /**
         * @brief Build a set of unique to unique paths
         * @param pathInfo Path information (input)
         * @param MCD Maximum conflict degree
         * @param un Unique nodes (input/output) -- can be modified (!)
         * @param u2u Unique-to-unique paths (output)
         */
        void buildUnique2Unique(const BiNodeMap<PathInfo>& pathInfo, float MCD,
                                const BiNodeMap<bool>& un,
                                NodeChainSet& u2u) const;

        /**
         * Filter conflicting and/or overlapping reductions
         * @param reductions Set of unique-to-unique reductions (input / output)
         * @param pathInfo Path information (input)
         * @param MCD Maximum conflict degree
         * @param relThreshold Filter relThreshold least frequent reductions
         */
        void filterReductions(NodeChainSet& red,
                              const BiNodeMap<PathInfo>& pathInfo,
                              float MCD, float relThreshold) const;

        /**
         * @brief Cluster u2u's into linear (or circular) contigs
         * @param u2u Unique-to-unique paths (inputs)
         * @param contigs Contigs (output)
         * @param node2contig Mapping between node and contig
         */
        void buildContigs(const NodeChainSet& u2u,
                          BiNodeMap<Contig>& contigs,
                          BiNodeMap<ContigID>& node2contig) const;

        std::vector<std::pair<ContigID, int>>
                getContigDest(NodeID srcID,
                              const BiNodeMap<PathInfo>& pathInfo,
                              BiNodeMap<bool>& un,
                              const BiNodeMap<Contig>& contigs,
                              const BiNodeMap<int>& node2contig) const;

        NodeChain connectContigs(ContigID srcContID, ContigID dstContID,
                                 const BiNodeMap<Contig>& contigs,
                                 const BiNodeMap<PathInfo>& pathInfo,
                                 float MCD) const;

        /**
         * Get the reductions from aligned reads
         * @param library Read library
         * @param nodeMap Multiplicity estimations for nodes
         * @param edgeMap Multiplicity estimations for edges
         * @param MCD Maximum conflict degree
         * @param minCount Minimum count (used to filter paired-end reads)
         * @param reductions Derived reductions (output)
         * @param looseEnds Loose ends
         */
        void getReductions(const LibraryContainer& library,
                           const NodeMap<int>& nodeMap,
                           const EdgeMap<int>& edgeMap,
                           float MCD, int minCount,
                           std::vector<NodeChain>& reductions,
                           std::vector<NodeChain>& looseEnds) const;

        /**
         * Extract the contigs from the graph, taking into account read info
         * @param library Read library (input)
         * @param nodeMap Multiplicity estimations for nodes (input)
         * @param edgeMap Multiplicity estimations for edges (input)
         * @param contigs Contigs (output)
         * @param contigName Contig names (output)
         */
        void getContigs(const LibraryContainer& library,
                        NodeMap<int>& nodeMult, EdgeMap<int>& edgeMap,
                        std::vector<NodeChain>& contigs,
                        std::vector<std::string>& contigName);

        /**
         * Get singleton contigs (one contig per graph node)
         * @param contigs Contigs (output)
         * @param contigName Contig names (output)
         */
        void getSingletonContigs(NodeMap<Multiplicity>& nodeMult,
                                 std::vector<NodeChain>& contigs,
                                 std::vector<std::string>& contigName) const;

        /**
         * Get the reductions from aligned paired-end reads
         * @param library Read library
         * @param nodeMap Multiplicity estimations for nodes
         * @param edgeMap Multiplicity estimations for edges
         * @param reductions Derived reductions (output)
         */
        void getPERReductions(const LibraryContainer& library,
                              NodeMap<Multiplicity>& nodeMap,
                              EdgeMap<Multiplicity>& edgeMap,
                              std::vector<NodeChain>& reductions);

        /**
         * Correct the node/edge multiplicity given the reductions and loose
         * ends. If the estimated node/edge multiplicity is smaller than
         * dictated by the reductions, adjust it accordingly
         * @param reductions Reductions
         * @param looseEnds Loose ends
         * @param nodeMult Estimated node multiplicity (input/output)
         * @param edgeMult Estimated edge multiplicity (input/output)
         * @param redNodeMult Node frequency in the reductions
         * @param redEdgeMult Edge frequency in the reductions
         */
        void adjustMult(const std::vector<NodeChain>& reductions,
                        const std::vector<NodeChain>& looseEnds,
                        NodeMap<int>& nodeMult, EdgeMap<int>& edgeMult,
                        NodeMap<int>& redNodeMult, EdgeMap<int>& redEdgeMult);

        /**
         * Apply a set of reductions to a de Bruijn graph
         * @param reductions Vector of reductions to apply
         * @param nodeMap Multiplicity of nodes in reductions (will be adjusted)
         * @param edgeMap Multiplicity of edges in reductions (will be adjusted)
         */
        void applyReductions(const std::vector<NodeChain>& reductions,
                             NodeMap<int>& nodeMap, EdgeMap<int>& edgeMap);

        /**
         * Get a list of nodes reachable from srcID as well as their shortest
         * pathlength from srcID (pathlength excludes both srcID and target)
         * @param srcID Source node identifier
         * @param nodeDist Map with <NodeID, pathlength> (output)
         * @param maxLen Maximum length (input)
         * @param maxNode Maximum number of nodes visited (input)
         */
        void getSubgraph(NodeID srcID, int maxLen, int maxNodes,
                         std::map<NodeID, int>& nodeDist);

        /**
         * Get a list of nodes between tangle srcID and dstID
         * @param srcID Source node identifier
         * @param dstID Destination node identifier
         * @param nodes Set with nodes in the subgraph (including srcID/dstID)
         */
        void getSubgraph(NodeID srcID, NodeID dstID, std::set<NodeID>& nodes);

        /**
         * Check if a (sub)graph has a unique Eulerian path (BEST theorem)
         * @param srcID Source node identifier
         * @param dstID Destination node identifier
         * @param edgeMult Edge multiplicities (only those >= 1 !)
         * @return True of false
         */
        bool hasUniqueEulerianPath(NodeID srcID, NodeID dstID,
                                   const std::map<EdgeID, int>& edgeMult) const;

        /**
         * Get Eulerian path from srcID to dstID
         * @param srcID Source node identifier
         * @param edgeMap Multiplicity of edges
         * @return Path from srcID to dstID
         */
        std::vector<NodeID> getEulerianPath(NodeID srcID,
                const std::map<EdgeID, int>& edgeMult) const;

        /**
         * Do all paths with maxLen that originate from srcID lead to dstID?
         * ("do all paths lead to Rome?")
         * @param srcID Source node identifier
         * @param dstID Destination node identifier
         * @param maxLen Maximum length
         * @return true of false
         */
        bool toRome(NodeID srcID, NodeID dstID, int maxLen);

        /**
         * Is there a path from srcID to dstID with avg. cov > specified
         * @param srcID Source node identifier
         * @param dstID Destination node identifier
         * @param minCov Minimum coverage
         * @param maxLen Maximum length
         * @return true of false
         */
        bool findPath(NodeID srcID, NodeID dstID, double minCov, int maxLen) const;

        /**
         * Find a unique path from srcID to dstID
         * @param srcID Source node identifier
         * @param dstID Destination node identifier
         * @param pf Path finder
         * @param MCD Maximum Confict Degree
         * @param path Output path
         * @return true of false
         */
        bool findPath(NodeID srcID, NodeID dstID,
                      const BiNodeMap<PathFinder> &pf,
                      float MCD, NodeChain& path) const;

        /**
         * Find a parallel path to a given tip (see definition for details)
         * @param tipID Tip identifier
         * @param maxLen Maximum length of the tip
         * @return Parallel path struct
         */
        ParallelPath findSimilarParallelTip(NodeID tipID, NodeLength maxLen);

        /**
         * Find a parallel path to a given node (see definition for details)
         * @param id Node identifier
         * @return Parallel path struct
         */
        ParallelPath findSimilarParallelPath(NodeID id);

        /**
         * Given a nodeID, find the shortest path that connects a unique
         * predecessor of nodeID with a unique ancestor of nodeID such that
         * each internal node of the path has higher coverage than nodeID
         * @param id Node identifier
         * @return Resulting path
         */
        std::vector<NodeID> findParallelPath(NodeID id);

        /**
         * Find tangles in the graph. These are pairs <srcID, dstID> such that
         * all paths from srcID *have* to exit through dstID and vice versa
         * @param nodeMult Node multiplicity estimates
         * @param tangles Tangles (output)
         */
        void findTangles(NodeMap<int>& nodeMult,
                         std::vector<NodePair>& tangles);

        /**
         * Try and resolve a tangle
         * @param tangle Tangle <srcID, dstID>
         * @param pathInfo Path information
         * @param nodeMult Node multiplicity estimates
         * @param edgeMult Edge multiplicity estimates
         * @param MCD Maximum conflict degree
         * @return Resolved tangle or empty if not successful
         */
        NodeChain resolveTangle(const NodePair& tangle,
                                BiNodeMap<PathInfo>& pathInfo,
                                NodeMap<int>& nodeMult, EdgeMap<int>& edgeMult,
                                double MCD);

        /**
         * Try and resolve the tangles
         * @param library Read library
         * @param nodeMult Node multiplicity estimates
         * @param edgeMult Edge multiplicity estimates
         * @param tangles Tangles
         * @param reductions Reductions (output)
         * @param redNodeMult Node multiplicity of reductions (output)
         * @param redEdgeMult Edge multiplicity of reductions (output)
         */
        void resolveTangles(LibraryContainer& library,
                            NodeMap<int>& nodeMult, EdgeMap<int>& edgeMult,
                            std::vector<NodePair>& tangles,
                            std::vector<NodeChain>& reductions,
                            NodeMap<int>& redNodeMult,
                            EdgeMap<int>& redEdgeMult);

        /**
         * Find linear paths of unique connected nodes
         * @param nodeMult Node multiplicity estimates
         * @param edgeMult Edge multiplicity estimates
         */
        void findLinearUniquePaths(NodeMap<int>& nodeMult,
                                   EdgeMap<int>& edgeMult);

        /**
         * Get node coverage histogram
         * @param hist < coverage, weight > histogram (output)
         */
        void getNodeCovHist(std::map<int, double>& hist);

        /**
         * Defragment graph
         */
        void defrag() {
                defragNodes();
                defragArcs();
        }

        /**
         * Reverse the orientation of a read alignment
         * @param nc Node chains
         * @param begR Begin position in read of alignment
         * @param endR End position in read of alignment
         * @param begAln Begin position in first node
         * @param endAln Begin position in last node
         * @param RL Read length
         */
        void reverseRead(NodeChain& nc, int& begR, int& endR,
                         int& begAln, int& endAln, int RL) const {
                if (nc.empty())
                        return;
                begR = RL - begR;
                endR = RL - endR;
                std::swap(begR, endR);
                begAln = getSSNode(nc.front()).length() - begAln;
                endAln = getSSNode(nc.back()).length() - endAln;
                std::swap(begAln, endAln);
                nc.revCompl();
        }

        /**
         * Get concatenation (if any) centered around seed node
         * @param seedID Node identifier of the seed node
         * @return Concatenation in the form a reduction
         */
        NodeChain getConcatenation(NodeID seedID);

        /**
         * Get reductions in the form of linear path concatenations
         * @param reductions Derived reductions (output)
         **/
        void getConcatenations(std::vector<NodeChain>& reductions);
};

#endif
