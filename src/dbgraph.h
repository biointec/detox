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
#include "kmernpp.h"

#include <mutex>
#include <vector>
#include <queue>

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class Settings;
class LibraryContainer;
class ReadRecord;
class Multiplicity;

// ============================================================================
// DIJKSTRA AUXILIARY CLASSES
// ============================================================================

class NodeDFS {

public:
        NodeID nodeID;          // current node identifier
        int depth;             // length to current node

        /**
         * Default constructor
         * @param nodeID node identifier
         * @param depth length to the current node
         */
        NodeDFS(NodeID nodeID, int depth) :
                nodeID(nodeID), depth(depth) {};
};

struct NodeDFSComp {

        /**
         * Compare by depth (use greater because priority_queue.top() returns
         * by default the greatest element. We want to return the smallest.)
         */
        bool operator()(const NodeDFS& f, const NodeDFS& s) {
                return f.depth > s.depth;
        }
};

class NodeRepDepth {

public:
        NodeRep nodeRep;        // current node identifier
        int depth;              // depth of the current node

        /**
         * Default constructor
         * @param nodeRep node identifier
         * @param depth depth of the current node
         */
        NodeRepDepth(NodeRep nodeRep, int depth) :
                nodeRep(nodeRep), depth(depth) {};
};

struct NodeRepComp {

        /**
         * Compare by depth (use greater because priority_queue.top() returns
         * by default the greatest element. We want to return the smallest.)
         */
        bool operator()(const NodeRepDepth& f, const NodeRepDepth& s) {
                return f.depth > s.depth;
        }
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
         **/
        void removeNode(NodeID nodeID);

        /**
         * Remove an arc between two nodes
         * @param leftID Identifier of the left node
         * @param rightID Identifier of the right node
         */
        void removeArc(NodeID leftID, NodeID rightID);

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
        void covCount(const ReadRecord& rr,
                      const KmerNPPTable& table);

        /**
         * Thread entry point for counting node and arc coverage
         * @param inputs Library container with input read files
         * @param table < Kmer, NPP > table with all contig endpoint kmers
         */
        void covCountThread(LibraryContainer& inputs,
                            const KmerNPPTable& table);

public:
        /**
         * Default constructor
         * @param settings Const-ref to settings object
         */
        DBGraph(const Settings& settings) : settings(settings), nodes(NULL),
                arcs(NULL), numNodes(0), numArcs(0) {};

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
         * Check whether a node exists
         * @param nr Node representative
         * @return true of false
         */
        bool nodeExists(NodeRep nr) const {
                if (nr.getNodeID() == 0)
                        return false;
                if (abs(nr.getNodeID()) > getNumNodes())
                        return false;
                return getSSNode(nr.getNodeID()).isValid();
        }

        /**
         * Check whether a edge exists
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

        /**
         * Check whether two NodePosPairs are consecutive in the graph
         * @param left Left NodePosPair
         * @param right Right NodePosPair
         * @return True if the NodePosPairs are consecutive
         */
        bool consecutiveNPP(const NodePosPair& left,
                            const NodePosPair& right) const;

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
                                 std::vector<std::pair<NodeID, NodeID> > edges,
                                 const std::map<NodeRep, Multiplicity>& estNodeMult,
                                 const std::map<EdgeRep, Multiplicity>& estEdgeMult,
                                 const std::map<NodeRep, int>& trueNodeMult,
                                 const std::map<EdgeRep, int>& trueEdgeMult) const;

        /**
         * Check the validity of the de Bruijn graph
         */
        void sanityCheck();

        /**
         * Get a subgraph of the de Bruijn graph
         * @param todo Priority queue with nodes to handle
         * @param nodes Set of nodes representatives in the subgraph (output)
         * @param edges Set of edge representatives in the subgraph (output)
         * @param maxDepth Maximum depth (in terms of number of nodes)
         */
        void getSubgraph(std::priority_queue<NodeRepDepth, std::vector<NodeRepDepth>,
                         NodeRepComp>& todo, std::set<NodeRep>& nodes,
                         std::set<EdgeRep>& edges, size_t maxDepth) const;

        /**
         * Get a subgraph of the de Bruijn graph
         * @param seedNode Seed node representative
         * @param nodes Set of nodes representatives in the subgraph (output)
         * @param edges Set of edge representatives in the subgraph (output)
         * @param maxDepth Maximum depth (in terms of number of nodes)
         */
        void getSubgraph(NodeRep seedNode, std::set<NodeRep>& nodes,
                         std::set<EdgeRep>& edges, size_t maxDepth = 0) const;

        /**
         * Get a subgraph of the de Bruijn graph
         * @param seedEdge Seed edge representative
         * @param nodes Set of nodes representatives in the subgraph (output)
         * @param edges Set of edge representatives in the subgraph (output)
         * @param maxDepth Maximum depth (in terms of number of nodes)
         */
        void getSubgraph(EdgeRep seedEdge, std::set<NodeRep>& nodes,
                         std::set<EdgeRep>& edges, size_t maxDepth = 0) const;

        /**
         * Get a directed subgraph of the de Bruijn graph
         * @param seedNode Seed node
         * @param nodes Vector of nodes in the subgraph (output)
         * @param edges Vector of edges in the subgraph (output)
         * @param maxDepth Maximum depth (in terms of number of nodes)
         */
        void getSubgraph(NodeID seedNode, std::vector<NodeID>& nodes,
                         std::vector<std::pair<NodeID, NodeID> >& edges,
                         size_t maxDepth = 0) const;

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
         * Get a nodes with an average coverage below the threshold
         * @param threshold Coverage threshold
         * @return A vector of node representations
         */
        std::vector<NodeRep> getLowCovNodes(double threshold) const;

        /**
         * Get tips with an average coverage below the threshold
         * @param threshold Coverage threshold
         * @return A vector of node representations
         */
        std::vector<NodeRep> getLowCovTips(double threshold) const;

        /**
         * Get bubbles with an average coverage below the threshold
         * @param threshold Coverage threshold
         * @return A vector of node representations
         */
        std::vector<NodeRep> getLowCovBubbles(double threshold) const;

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
         * Convert a vector of overlapping nodes to a string
         * @param nodeSeq A deque of overlapping nodes
         * @param output An stl string (output)
         */
        void convertNodesToString(const std::vector<NodeID> &nodeSeq,
                                  std::string &output);

        /**
         * Concatentate linear paths
         * @return True if at least one node was merged
         **/
        bool concatenateNodes();

        /**
         * Remove of list of nodes from the graph
         * @param nodes Vector of node representatives to remove
         */
        void removeNodes(const std::vector<NodeRep>& nodes);

        /**
         * Remove of list of edges from the graph
         * @param nodes Vector of edge representatives to remove
         */
        void removeEdges(const std::vector<EdgeRep>& edges);

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
};

#endif
