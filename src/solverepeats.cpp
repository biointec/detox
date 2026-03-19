/******************************************************************************
 *   Copyright (C) 2014 - 2022 Jan Fostier (jan.fostier@ugent.be)             *
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

#include "dbgraph.h"
#include "pathinfo.h"
#include "crfmult.h"

using namespace std;

void DBGraph::addReadInfo(const NodeChain& read,
                          const BiNodeMap<bool>& noi,
                          BiNodeMap<ChainColl>& chainColl) const
{
        for (size_t i = 0; i < read.size(); i++) {
                if (!get(noi, read[i]))
                        continue;       // not unique
                NodeChain readSuffix(read.begin() + i, read.end(), 1);
                if (readSuffix.size() <= 1)
                        continue;       // uninformative
                chainColl[read[i]].addChain(readSuffix);
        }
}

void DBGraph::addReadInfo(const LibraryContainer& library,
                          const BiNodeMap<bool>& noi,
                          BiNodeMap<PathInfo>& readInfo) const
{
        BiNodeMap<ChainColl> chainColl;

        // process all read traces
        for (const Library& lib : library) {
                const auto [fn1, fn2] = lib.getAlnFilename();

                ifstream ifs1(fn1);
                ifstream ifs2(fn2);
                string line1, line2;
                while (getline(ifs1, line1) && getline(ifs2, line2)) {
                        NodeChain L, R;

                        istringstream iss(line1);
                        NodeID nodeID;
                        size_t ignore;
                        iss >> ignore >> ignore;
                        while (iss >> nodeID)
                                L.push_back(nodeID);
                        L.setCount(1);

                        iss = istringstream(line2);
                        iss >> ignore >> ignore;
                        while (iss >> nodeID)
                                R.push_back(nodeID);
                        R.setCount(1);

                        // build the path info for L and R
                        addReadInfo(L, noi, chainColl);
                        addReadInfo(R, noi, chainColl);

                        // do the same for the reverse complement of L and R
                        L.revCompl();
                        R.revCompl();
                        addReadInfo(L, noi, chainColl);
                        addReadInfo(R, noi, chainColl);
                }       // of paired-end reads
        }       // of libraries

        // finish adding paths
        readInfo.clear();
        for (auto& [nodeID, cc] : chainColl) {
                readInfo[nodeID] = PathInfo(cc);
                cc.clear();
        }
}

void DBGraph::addPERInfo(const NodeChain& L, const NodeChain& R,
                         size_t readID, const BiNodeMap<bool>& noi,
                         BiNodeMap<PathInfo>& readInfo) const
{
        set<NodeID> toAdd;
        for (size_t i = 0; i < L.size(); i++)
                if (get(noi, L[i]))
                        toAdd.insert(L[i]);

        for (size_t i = 0; i < R.size(); i++)
                if (get(noi, R[i]))
                        toAdd.insert(R[i]);

        for (NodeID nodeID : toAdd) {
                size_t i = 0;
                for ( ; i < L.size(); i++)
                        if (L[i] == nodeID)
                                break;

                size_t j = 0;
                for ( ; j < R.size(); j++)
                        if (R[j] == nodeID)
                                break;
                if (j == R.size())
                        j = 0;  // if R does not contain noi, add R completely

                NodeChain LSuffix(L.begin() + i, L.end(), 1);
                NodeChain RSuffix(R.begin() + j, R.end(), 1);

                readInfo[nodeID].addPer(LSuffix, RSuffix, readID);
        }
}

void DBGraph::addPERInfo(const LibraryContainer& library,
                         const BiNodeMap<bool>& noi,
                         BiNodeMap<PathInfo>& readInfo) const
{
        // process all read traces
        size_t readID = 0;      // unique read identifier
        for (const Library& lib : library) {
                const auto [fn1, fn2] = lib.getAlnFilename();

                ifstream ifs1(fn1);
                ifstream ifs2(fn2);
                string line1, line2;
                while (getline(ifs1, line1) && getline(ifs2, line2)) {
                        NodeChain L, R;

                        istringstream iss(line1);
                        NodeID nodeID;
                        size_t ignore;
                        iss >> ignore >> ignore;
                        while (iss >> nodeID)
                                L.push_back(nodeID);
                        L.setCount(1);

                        iss = istringstream(line2);
                        iss >> ignore >> ignore;
                        while (iss >> nodeID)
                                R.push_back(nodeID);
                        R.setCount(1);

                        // build the path info for L and R
                        addPERInfo(L, R, readID, noi, readInfo);

                        // do the same for the reverse complement PER
                        swap(L, R);
                        L.revCompl();
                        R.revCompl();
                        addPERInfo(L, R, readID, noi, readInfo);

                        readID++;
                }
        }
}

void DBGraph::filterPathInfo(BiNodeMap<PathInfo>& pathInfo,
                             float MCD, int minCount) const
{
        for (auto& [ignore, pi] : pathInfo)
                pi.filterDirPaths(MCD);
}

void DBGraph::buildUnique2Unique(const BiNodeMap<PathInfo>& pathInfo,
                                 float MCD, const BiNodeMap<bool>& un,
                                 NodeChainSet& u2u) const
{
        u2u.clear();

        BiNodeMap<ConsNodeChain> MLP;
        for (const auto& [srcID, srcPI] : pathInfo)
                MLP[srcID] = srcPI.getMostLikelyPath(MCD);

        for (const auto& [srcID, mlp] : MLP) {
                for (size_t i = 1; i < mlp.size(); i++) {
                        NodeID dstID = mlp[i].first;
                        if (!get(un, -dstID))
                                continue;       // not unique
                        NodeChain nc = ConsNodeChain(mlp.begin(),
                                mlp.begin() + i + 1).toNodeChain();
                        NodeChain rc = nc.getRevCompl();

                        if (MLP.find(-dstID) == MLP.end())
                                continue;
                        if (!MLP[-dstID].isPrefixed(rc))
                                continue;
                        u2u.insert(nc);
                        u2u.insert(rc);
                }
        }
}

void DBGraph::filterReductions(NodeChainSet& red,
                               const BiNodeMap<PathInfo>& readInfo,
                               float MCD, float relThreshold) const
{
        NodeChainSet u2u;
        swap(u2u, red);

        // flag all unique nodes
        NodeMap<bool> un;
        for (const NodeChain& nc : u2u) {
                un[nc.front()] = true;
                un[-nc.back()] = true;
        }

        // flag nodes with conflicting u2u's
        BiNodeMap<ConsNodeChain> nodePath;
        for (const NodeChain& nc : u2u) {
                for (size_t i = 0; i < nc.size(); i++) {
                        NodeID id = nc[i];
                        if (!get(un, id))
                                continue;
                        NodeChain suffix(nc.begin() + i, nc.end());

                        // Important to check the suffix against the readInfo!
                        // In case mult(id) > 1, no u2u might originate from
                        // id and then the nodePath check below will not fire!
                        if (!get(readInfo, id).isCompatible(suffix, MCD))
                                un[id] = un[-id] = false;       // conservative

                        // Check against the other u2u's
                        if (!nodePath[id].mergeNodeChain(suffix))
                                un[id] = un[-id] = false;       // conservative
                }
        }

        // remove paths between nodes that are no longer unique
        for (auto it = u2u.begin(); it != u2u.end(); ) {
                if (!un[it->front()] || !un[-it->back()])
                        it = u2u.erase(it);
                else
                        ++it;
        }

        // compute a confidence threshold
        map<size_t, size_t> u2uHistogram;
        for (const NodeChain& nc : u2u)
                u2uHistogram[nc.getCount()]++;

        size_t thresHold = 0, cumObs = 0;
        for (const auto [count, numObs] : u2uHistogram) {
                cumObs += numObs;
                if ((float)cumObs / (float)u2u.size() < relThreshold)
                        thresHold = count;
                else
                        break;
        }

        thresHold = max<size_t>(1, thresHold);

        // build a set of shortest reductions
        red.clear();
        size_t numFiltered = 0;
        for (const NodeChain& nc : u2u) {
                if (nc.getCount() <= thresHold) {
                        numFiltered++;
                        continue;
                }
                for (size_t b = 0, i = 1; i < nc.size(); i++) {
                        if (!get(un, nc[i]))
                                continue;
                        NodeChain r(nc.begin() + b, nc.begin() + i + 1,
                                    nc.getCount());
                        b = i;

                        red.insert(r);
                        r.revCompl();
                        red.insert(r);
                }
        }

        cout << "Filtered " << numFiltered << " reductions that have fewer "
             << "than " << thresHold + 1 << " reads supporting them." << endl;
}

void DBGraph::buildContigs(const NodeChainSet& red,
                           BiNodeMap<Contig>& contigs,
                           BiNodeMap<ContigID>& node2contig) const
{
        contigs.clear();
        node2contig.clear();    // initially: all nodes are unassigned
        for (const auto& nc : red)
                node2contig[nc.front()] = node2contig[-nc.front()] = 0;

        int contigID = 1;
        for (const NodeChain& nc : red) {
                if (node2contig[nc.front()] != 0)
                        continue;       // already handled, skip

                // create a new contig
                ConsNodeChain& contig = contigs[contigID];
                node2contig[ nc.front()] = node2contig[ nc.back()] =  contigID;
                node2contig[-nc.front()] = node2contig[-nc.back()] = -contigID;
                contig = ConsNodeChain(nc);

                for (int i = 0; i < 2; i++) {   // extend in both directions
                        while (true) {
                                auto it = red.begin(contig.back().first);
                                if (it == red.end(contig.back().first))
                                        break;  // cannot extend further
                                const NodeChain& ext = *it;
                                const int c = ext.getCount();
                                transform(ext.begin() + 1, ext.end(),
                                          std::back_inserter(contig),
                                          [c](const NodeID& id){
                                                return make_pair(id, c); });
                                if (contig.front().first == contig.back().first)
                                        break;  // circular contig reached
                                node2contig[ ext.back()] =  contigID;
                                node2contig[-ext.back()] = -contigID;
                        }
                        if (contig.front().first == contig.back().first)
                                break;  // circular contig reached
                        // swap direction of contig
                        contig.revCompl();
                        contigID = -contigID;
                }

                contigs[-contigID] = contigs[contigID].getRevCompl();
                contigID++;
        }
}

vector<pair<ContigID, int>>
DBGraph::getContigDest(ContigID srcContigID, const BiNodeMap<PathInfo>& pathInfo,
                       BiNodeMap<bool>& un, const BiNodeMap<Contig>& contigs,
                       const BiNodeMap<int>& node2contig) const
{
        const Contig& srcCont = get(contigs, srcContigID);

        // count how many nodes at the end of srcContig are taken into account
        size_t len = 0, bIdx = srcCont.size();
        while (bIdx > 0) {
                bIdx--;
                len += getSSNode(srcCont[bIdx].first).getMarginalLength();
                if (len > MAX_CONTIG_ANCHOR_LEN)
                        break;
        }

        map<ContigID, int> dstContigID;                 // <dstContigID, count>
        map<ContigID, ReadIDSet> dstContigReadIDs;      // <dstContigID, readIDs>

        for (size_t i = bIdx; i < srcCont.size(); i++) {
                NodeID srcID = srcCont[i].first;
                if (!un[srcID])
                        continue;
                const PathInfo& pi = get(pathInfo, srcID);
                BiNodeMap<set<size_t>> dstPer = pi.getDstReadSet(un);

                for (const auto& [dstID, readIDs] : dstPer) {
                        if (abs(get(node2contig, dstID)) == abs(srcContigID))
                                continue;

                        ContigID thisDstContID = get(node2contig, dstID);
                        for (auto rID : readIDs) {
                                if (dstContigReadIDs[thisDstContID].insert(rID).second)
                                        dstContigID[thisDstContID]++;
                        }
                }
        }

        auto sortByWeight = [] (pair<ContigID, int> const& s1,
                                pair<ContigID, int> const& s2) -> bool {
                return s1.second > s2.second;
        };

        vector<pair<ContigID, int>> result;
        copy(dstContigID.begin(), dstContigID.end(), back_inserter(result));
        sort(result.begin(), result.end(), sortByWeight);

        return (result.size() <= 5) ? result:
                vector<pair<ContigID, int>>(result.begin(), result.begin() + 5);

}

bool DBGraph::findPath(NodeID srcID, NodeID dstID,
                       const BiNodeMap<PathFinder> &pf,
                       float MCD, NodeChain& path) const
{
        // sanity checks
        assert(NodeRep(srcID) != NodeRep(dstID));
        assert(nodeExists(srcID));
        assert(nodeExists(dstID));

        path.clear();
        path.setCount(numeric_limits<int>::max());
        path.push_back(srcID);

        // find a non-conflicting path from srcPath.back() to dstPath[dstIdx]
        while (true) {
                NodeID currID = path.back();
                SSNode n = getSSNode(currID);

                NodeID nextID = 0;
                if (n.numRightArcs() == 0)
                        return false;   // no right arc
                else if (n.numRightArcs() == 1)
                        nextID = n.rightBegin()->getNodeID();
                else
                        nextID = get(pf, currID).getNextID(path, MCD);

                if (nextID == 0) {      // process the right arcs
                        Coverage maxVal = 0, sumVal = 0;
                        for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                                sumVal += it->getCov();
                                if (it->getCov() > maxVal) {
                                        nextID = it->getNodeID();
                                        maxVal = it->getCov();
                                }
                        }

                        if ((float(sumVal - maxVal) / float(sumVal)) > MCD)
                                return false;   // conflict degree exceeded
                }

                if (path.size() > 500)
                        return false;   // stuck in a loop

                path.setCount(min<int>(path.getCount(), getArc(currID, nextID).getCov()));
                path.push_back(nextID);

                if (path.back() == dstID)
                        return true;
        }
}

NodeChain DBGraph::connectContigs(ContigID srcContID, ContigID dstContID,
                                  const BiNodeMap<Contig>& contigs,
                                  const BiNodeMap<PathInfo>& readInfo,
                                  float MCD) const
{
        // get source and destination node
        NodeID srcID = get(contigs, srcContID).back().first;
        NodeID dstID = get(contigs, dstContID).front().first;

        // get a set of reads that pass through either src or dst nodes
        vector<NodeChain> reads; set<size_t> readIDs;
        get(readInfo,  srcID).addUniqueReads(reads, readIDs, false);
        get(readInfo, -dstID).addUniqueReads(reads, readIDs, true);
        readIDs.clear();        // readIDs avoids duplicates, no longer needed

        // flatten in a list of nodes and edges
        NodeMap<int> nodes; EdgeMap<int> edges;
        for (const auto& r : reads) {
                for (size_t i = 0; i < r.size(); i++) {
                        nodes[NodeRep(r[i])]++;
                        if (i > 0)
                                edges[EdgeRep(r[i-1], r[i])]++;
                }
        }

        // build a bi-directed de Bruijn subgraph
        BiNodeMap<NodeID> subToOrig, origToSub;
        DBGraph subGraph(*this, nodes, edges, subToOrig, origToSub);

        // figure out the node ID of srcID in the subgraph
        NodeID subSrcID = 0, subDstID = 0, counter = 0;
        for (auto [nr, ignore] : nodes) {
                counter++;
                if (nr == NodeRep(srcID))
                        subSrcID = (srcID > 0) ? counter : -counter;
                if (nr == NodeRep(dstID))
                        subDstID = (dstID > 0) ? counter : -counter;
        }

        BiNodeMap<NodePosPair> origToSub2;
        createNodeMapping(subGraph, srcID, subSrcID, origToSub2);

        vector<NodeChain> readSub;
        for (const auto& nc : reads)
                readSub.emplace_back(convertNC(nc, subGraph, origToSub2));

        BiNodeMap<PathFinder> pf;
        vector<NodeID> bn = subGraph.getBranchingNodes();
        for (auto e : bn)
                pf[e] = PathFinder();

        for (const auto& nc : readSub) {
                for (size_t i = 1; i + 1 < nc.size(); i++) {
                        if (pf.find(nc[i]) == pf.end())
                                continue;       // not of interest
                        pf[nc[i]].addChain(NodeChain(nc.begin(), nc.begin() + i), nc[i+1]);
                }
        }

        for (auto& [n, p] : pf)
                p.buildIndex();

        NodeChain nc; ConsNodeChain subNc;
        if (subGraph.findPath(origToSub[srcID], origToSub[dstID], pf, MCD, nc)) {
                ConsNodeChain temp(nc);
                subNc = subGraph.convertCNC(temp, *this, subToOrig);
                //cout << "Found path on the uncleaned subgraph: " << subNc << endl;
                return subNc.toNodeChain();
        }

        // prevent src and dst node from begin removed / concatenated
        set<NodeRep> toKeep = {origToSub[srcID], origToSub[dstID]};

        // graph cleaning
        subGraph.concatenateNodes(toKeep);
        while (subGraph.removeLowCovTips(numeric_limits<double>::max(),
                                         MAX_TIP_LEN, toKeep, false) > 0);
        while (subGraph.removeLowCovBubbles(numeric_limits<double>::max(),
                                            MAX_BUBBLE_LEN, toKeep) > 0);

        origToSub2.clear();
        createNodeMapping(subGraph, srcID, subSrcID, origToSub2);

        readSub.clear();
        for (const auto& nc : reads)
                readSub.emplace_back(convertNC(nc, subGraph, origToSub2));

        pf.clear();
        bn = subGraph.getBranchingNodes();
        for (auto e : bn)
                pf[e] = PathFinder();

        for (const auto& nc : readSub) {
                for (size_t i = 1; i + 1 < nc.size(); i++) {
                        if (pf.find(nc[i]) == pf.end())
                                continue;       // not of interest
                        pf[nc[i]].addChain(NodeChain(nc.begin(), nc.begin() + i), nc[i+1]);
                }
        }

        for (auto& [n, p] : pf)
                p.buildIndex();

        if (subGraph.findPath(origToSub[srcID], origToSub[dstID], pf, MCD, nc)) {
                ConsNodeChain temp(nc);
                subNc = subGraph.convertCNC(temp, *this, subToOrig);
                //cout << "Found path on the cleaned subgraph: " << subNc << endl;
                return subNc.toNodeChain();
        }

        // ============== CYTOSCAPE ================
        /*cout << "Writing Cytoscape graph..." << endl;
        vector<NodeID> cytNodes; vector<EdgeID> cytEdges;
        //subGraph.getGraph(cytNodes, cytEdges);
        //subGraph.writeCytoscapeGraph("graphPER", cytNodes, cytEdges, NodeMap<Multiplicity>(), EdgeMap<Multiplicity>(), NodeMap<int>(), EdgeMap<int>());

        subGraph.getGraph(cytNodes, cytEdges);
        subGraph.writeCytoscapeSubgraph("graphPER", cytNodes, cytEdges, subToOrig);
        exit(0);*/
        // ============== CYTOSCAPE ================

        return NodeChain();
}

void DBGraph::buildPerInfo(const LibraryContainer& library,
                           const BiNodeMap<Contig>& contigs,
                           const BiNodeMap<int>& node2contig,
                           BiNodeMap<PathInfo>& readInfo, float MCD,
                           NodeChainSet& u2uPer) const
{
        // get all nodes that can be used to identify a contig
        BiNodeMap<bool> un;
        for (const auto& [id, ignore] : node2contig)
                un[id] = true;

        addPERInfo(library, un, readInfo);

        for (const auto& [srcContID, srcContig] : contigs)
        //auto srcContID = -679;
        //const auto& srcContig = get(contigs, srcContID);
        {
                float bestScore = 0;
                vector<pair<ContigID, int>> dstContID_v =
                        getContigDest(srcContID, readInfo, un, contigs, node2contig);

                //cout << "Handling " << srcContID << "...";

                NodeID srcID = get(contigs, srcContID).back().first;
                const PathInfo& srcPI    = get(readInfo,  srcID);

                for (const auto [dstContID, perCount] : dstContID_v) {

                        // This is an approximation. It is possible that
                        // the overal score exceeds perCount
                        if (perCount < MCD * bestScore)
                                break;

                        if (srcContID > -dstContID)
                                continue;       // reverse-complement handled

                        //cout << "Handling" << srcContID << " to " << dstContID << endl;
                        //Util::startChrono();

                        /*cout << "Connecting contigs " << srcContID << " and " << dstContID
                             << "\n\t srcID: " << srcID << " - dstID: " << get(contigs, dstContID).front().first
                             << " (w = " << perCount << ")\n";*/

                        NodeChain nc = connectContigs(srcContID, dstContID,
                                                      contigs, readInfo, MCD);

                        //cout << Util::stopChronoStr() << endl;
                        if (nc.empty())
                                continue;

                        //cout << "Found path: " << nc << endl;

                        int ncCount = nc.getCount();
                        float geoMean = sqrt(float(ncCount) * float(perCount));
                        nc.setCount(min<int>(perCount, int(geoMean + 0.5f)));
                        bestScore = nc.getCount();
                        u2uPer.insert(nc);
                        nc.revCompl();
                        u2uPer.insert(nc);
                }
        }
}

void DBGraph::getReductions(const LibraryContainer& library,
                            const NodeMap<int>& nm,
                            const EdgeMap<int>& em,
                            float MCD, int minCount,
                            vector<NodeChain>& red,
                            vector<NodeChain>& looseEnds) const
{
        BiNodeMap<bool> un;             // indicate the unique nodes
        for (const auto& [id, m] : nm) {
                un[id] =  (m == 1);
                un[-id] = (m == 1);
        }

        // build the path info from the aligned reads
        BiNodeMap<PathInfo> pathInfo;
        addReadInfo(library, un, pathInfo);
        filterPathInfo(pathInfo, MCD, minCount);

        // filter nodes that loop onto themselves
        for (auto& [nodeID, pi] : pathInfo) {
                NodeChain nc = pi.getMostLikelyPath(MCD).toNodeChain();
                for (size_t i = 1; i < nc.size(); i++) {
                        if (NodeRep(nc[i]) == NodeRep(nodeID)) {
                                un[ nodeID] = false;
                                un[-nodeID] = false;
                                pi.clear();
                                pathInfo[-nodeID].clear();
                                break;
                        }
                }
        }

        // try and correct for nodes with estimated multiplicity zero that
        // may be part of a unique path anyhow
        /*for (auto& [thisID, pi] : pathInfo) {
               if (thisID != -3939 && thisID != 3939)
                        continue;
                ConsNodeChain cnc = pi.getMostLikelyPath(MCD);
                //cout << cnc << endl;
                for (auto [id, count] : cnc) {
                        if (get(nm, id) == 0 && count >= minCount) {
                                un[id]  = true;
                                un[-id] = true;
                                cout << "Setting " << id << " to unique" << endl;
                        }
                }
        }



        pathInfo.clear();
        addReadInfo(library, un, pathInfo);
        filterPathInfo(pathInfo, MCD, minCount);*/

        // build a set unique-to-unique reductions
        NodeChainSet u2u;
        Util::startChrono();
        buildUnique2Unique(pathInfo, MCD, un, u2u);

        // build a set of shortest (= minimal), non-conflicting reductions
        NodeChainSet minimalRed = u2u;
        filterReductions(minimalRed, pathInfo, MCD, 0.02f);
        cout << "Found " << minimalRed.size()/2 << " non-conflicting, "
             << "minimal reductions in " << Util::stopChronoStr() << "\n";

        // aggregate the reductions into larger contigs
        BiNodeMap<Contig> contigs;
        BiNodeMap<int> node2contig;

        buildContigs(minimalRed, contigs, node2contig);
        cout << "Clustered reductions into " << contigs.size()/2 << " contigs\n";

        NodeMap<bool> inContig;
        for (const auto& [ignore, c] : contigs)
                for (auto [id, ignore] : c)
                        inContig[id] = true;

        // also create contigs from isolated nodes
        int contigID = contigs.size() / 2 + 1;
        for (const auto [id, flag] : un) {
                if (!flag || id != NodeRep(id))
                        continue;
                if (!pathInfo[id].isUnique(MCD) || !pathInfo[-id].isUnique(MCD))
                        continue;       // extra check based on pathInfo
                if (inContig.find(id) != inContig.end())
                        continue;       // already in a contig

                node2contig[ id] =  contigID;
                node2contig[-id] = -contigID;
                contigs[ contigID] = ConsNodeChain(vector<NodeID>(1,  id));
                contigs[-contigID] = ConsNodeChain(vector<NodeID>(1, -id));
                contigID++;
        }

        /*cout << "Is unique: -7291: " << (un[-7291] == true ? "yes" : "no")  << endl;
        cout << "Is unique: 6335: " << (un[6335] == true ? "yes" : "no")  << endl;*/
        //cout << "CONTIG: " << node2contig[-9663] << "\t" << contigs[node2contig[-9663]] << endl;
        //cout << "CONTIG: " << node2contig[-2998] << "\t" << contigs[node2contig[-2998]] << endl;

        //exit(0);

        cout << "After adding isolated nodes: " << contigs.size()/2 << " contigs\n";

        // build a comprensive set of paired-end contig-to-contig connections
        NodeChainSet perLinks;
        buildPerInfo(library, contigs, node2contig, pathInfo, MCD, perLinks);

        /*for (auto& nc : u2u)
                if (nc.contains(3592) || nc.contains(3592))
                        cout << "U2U DIR: " << nc << endl;

        for (auto& nc : perLinks)
                if (nc.contains(3592) || nc.contains(3592))
                        cout << "U2U PER: " << nc << endl;*/

        for (const auto& nc : u2u)
                perLinks.insert(nc);

        BiNodeMap<ChainColl> cc;
        for (auto& nc : perLinks) {
                for (size_t i = 0; i < nc.size() - 1; i++) {
                        if (node2contig.find(nc[i]) != node2contig.end()) {
                                NodeChain suffix(nc.begin() + i, nc.end());
                                suffix.setCount(nc.getCount());
                                cc[nc[i]].addChain(suffix);
                        }
                }
        }

        BiNodeMap<PathInfo> perInfo;
        for (auto& [id, nc] : cc)
                perInfo[id] = PathInfo(nc);

        // build a set of unique-to-unique paired-end reductions
        u2u.clear();
        buildUnique2Unique(perInfo, MCD, un, u2u);

        // It is important to again bring in *all* u2us, not only the shortest,
        // filtered ones (minimalRed). Collectively with the u2uPers, this
        // could lead to new insights.
        minimalRed = u2u;

        // Pass in empty path info: we no longer filter on path info, since
        // findContig might have
        // used a simplified graph
        BiNodeMap<PathInfo> empty;
        filterReductions(minimalRed, empty, MCD, 0.02f);

        cout << "Number of minimal reductions (incl. paired-end): "
             << minimalRed.size()/2 << "\n";

        buildContigs(minimalRed, contigs, node2contig);
        cout << "Clustered reductions into " << contigs.size()/2 << " contigs\n";

        //cout << "CONTIG: " << node2contig[-9663] << "\t" << contigs[node2contig[-9663]] << endl;
        //cout << "CONTIG: " << node2contig[-4915] << "\t" << contigs[node2contig[-4915]] << endl;
        //cout << "CONTIG: " << node2contig[-2998] << "\t" << contigs[node2contig[-2998]] << endl;

        red.clear();
        for (auto& [ignore, cnc] : contigs) {
                NodeChain nc = cnc.toNodeChain();
                if (nc.size() <= 1)
                        continue;
                if (nc < nc.getRevCompl())
                        red.emplace_back(nc);
        }

        // label nodes that are part of a reduction
        BiNodeMap<bool> endPoint = un;
        for (const auto& nc : red) {
                assert(!nc.empty());    // sanity check
                endPoint[nc.front()] = false;
                for (size_t i = 1; i < nc.size() - 1; i++) {
                        endPoint[ nc[i]] = false;
                        endPoint[-nc[i]] = false;
                }
                endPoint[-nc.back()] = false;
        }

        // look at the unique nodes for which there is no reduction
        looseEnds.clear();
        for (const auto& [id, flag] : endPoint) {
                if (!flag)
                        continue;       // not a unique node or not an endpoint
                const PathInfo& pi = get(pathInfo, id);
                NodeChain le = pi.getMostLikelyPath(MCD).toNodeChain();
                for (size_t i = 1; i < le.size(); i++) {
                        if (!un[le[i]])
                                continue;
                        le = NodeChain(le.begin(), le.begin() + i);
                }

                if (le.size() > 1)
                        looseEnds.push_back(le);
        }

        cout << "Number of reductions: " << red.size() << endl;
        cout << "Number of loose ends: " << looseEnds.size() << endl;
}
