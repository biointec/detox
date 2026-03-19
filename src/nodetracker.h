/******************************************************************************
 *   Copyright (C) 2022 Jan Fostier (jan.fostier@ugent.be)                    *
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

#ifndef NODETRACKER_H
#define NODETRACKER_H

#include <vector>
#include "ssnode.h"

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class DBGraph;
class NodeChain;

// ============================================================================
// CLASS NODE INFO
// ============================================================================

class NodeInfo {

public:
        NodeID nodeID;          // node identifier
        NodeLength offset;      // node offset
        NodeLength margLen;     // node marginal length

        /**
         * Default constructor
         */
        NodeInfo(NodeID nodeID = 0, NodeLength offset = 0,
                 NodeLength margLen = 0) :
                nodeID(nodeID), offset(offset), margLen(margLen) {}
};

// ============================================================================
// CLASS NODE TRACKER
// ============================================================================

/*
 The Node Tracker keeps tracks of all nodes that will be deleted/modified
 because of graph reductions. This includes unique nodes that will become
 part of bigger nodes (and deleted) as well as repeated nodes that will
 be deleted. For unique nodes, we keep track of the offset within the new node
 as well as the original node length. For (deleted) repeat nodes, we keep
 track track of the original node length.

 The goal of this class is to provide functionality to update node alignments,
 taking into account
*/

class NodeTracker {

private:
        const DBGraph& dBG;             // const-ref to de Bruijn graph
        BiNodeMap<NodeInfo> log;        // node tracker

        /**
         * Log that node oldID has been removed and now exists within
         * node newID at a certain offset. Also log length of old node.
         * @param oldID Identifier of the old node
         * @param newID Identifier of the new node
         * @param offset Offset of oldID within newID
         * @param margLen Marginal node length of the old node
         */
        void logEvent(NodeID oldID, NodeID newID,
                      NodeLength offset, NodeLength margLen) {
                log[oldID] = NodeInfo(newID, offset, margLen);
        }

        /**
         * Get the new node identifier
         * @param oldID Old node identifier
         * @return New node identifier
         */
        NodeID getPresentNodeID(NodeID oldID) const {
                auto it = log.find(oldID);
                return (it == log.end()) ? oldID : it->second.nodeID;
        }

        /**
         * Get the node length of oldID (that may have been deleted)
         * @param oldID Old node identifier
         * @return Node length
         */
        NodeLength getNodeLength(NodeID oldID) const;

        /**
         * Get the node offset of oldID within newID
         * @param oldID Old node identifier
         * @return Node offset
         */
        NodeLength getNodeOffset(NodeID oldID) const {
                const auto it = log.find(oldID);
                return (it == log.end()) ? 0 : it->second.offset;
        }

        /**
         * An anchor node is a *unique* node that is now part of a bigger node
         * @param oldID Old node identifier
         * @return true or false
         */
        bool isAnchorNode(NodeID oldID) const {
                const auto it = log.find(oldID);
                return (it != log.end() && it->second.nodeID != 0);
        }

public:
        /**
         * Constructor
         * @param dBG Const-ref to de Bruijn graph
         */
        NodeTracker(const DBGraph& dBG) : dBG(dBG) {}

        /**
         * Node tracker empty()
         * @return true or false
         */
        bool empty() const {
                return log.empty();
        }

        /**
         * Build the node tracker (call BEFORE doing actual reductions)
         * @param reductions Graph reductions
         * @param nodeMult Estimated node multiplicity
         * @param redNodeMult Node frequency in the reductions
         */
        void build(const std::vector<NodeChain>& reductions,
                   const NodeMap<int>& nodeMult,
                   const NodeMap<int>& redNodeMult);

        /**
         * Build the node tracker without node multiplicity info
         * @param reductions Graph reductions
         */
        void build(const std::vector<NodeChain>& reductions);

        /**
         * Update a node chain taking into account graph manipulations.
         * (call AFTER doing actual reductions)
         * @param nc Node chain (input/output)
         * @param begAln Begin offset relative to the first node (input/output)
         * @param endAln End offset relative to the last node (input/output)
         */
        void updateNodeChain(NodeChain& nc, int& begAln, int& endAln) const;
};

#endif
