# /******************************************************************************
#  *   Copyright (C) 2014 - 2026 Aranka Steyaert (aranka.steyaert@ugent.be)     *
#  *   This file is part of Detox                                               *
#  *                                                                            *
#  *   This program is free software; you can redistribute it and/or modify     *
#  *   it under the terms of the GNU Affero General Public License as published *
#  *   by the Free Software Foundation; either version 3 of the License, or     *
#  *   (at your option) any later version.                                      *
#  *                                                                            *
#  *   This program is distributed in the hope that it will be useful,          *
#  *   but WITHOUT ANY WARRANTY; without even the implied warranty of           *
#  *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
#  *   GNU General Public License for more details.                             *
#  *                                                                            *
#  *   You should have received a copy of the GNU Affero General Public License *
#  *   along with this program; if not, see <https://www.gnu.org/licenses/>.    *
#  ******************************************************************************/
import random
import numpy as np
import copy
import datetime
import os
import argparse
import subprocess
import sys

import networkx as nx

import matplotlib.pyplot as plt
from collections import deque

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.nn import Linear

from torchinfo import summary
from torch.utils.tensorboard import SummaryWriter

from torch.utils.data import random_split

from sklearn.model_selection import StratifiedKFold

import seaborn as sns

# %%
# Dataset as used for the paper CHANGE THIS TO YOUR OWN DATA
data_folder = '../data2/'
directories = ['ECOLI',
               'EFAECALIS', 'EFAECALIS2',
               'MTUBERCULOSIS', 'MTUBERCULOSIS2', 'MTUBERCULOSIS3', 'MTUBERCULOSIS4', 'MTUBERCULOSIS5', 'MTUBERCULOSIS6',
               'NMENINGITIDIS',
               'PAERUGINOSA', 'PAERUGINOSA2', 'PAERUGINOSA3', 'PAERUGINOSA4', 'PAERUGINOSA5',
               'SAUREUS', 'SAUREUS2',
               'SENTERICA', 'SENTERICA2', 'SENTERICA3', 'SENTERICA4', 'SENTERICA5',
               'SHIGELLA']

groups = ['ECOLI',
          'EFAECALIS',
          'MTUBERCULOSIS',
          'NMENINGITIDIS',
          'PAERUGINOSA',
          'SAUREUS',
          'SENTERICA',
          'SHIGELLA']

training_subset = '/REAL/30/' # should contain truemult.node and truemult.edge files
# %%
class Graph:
    def __init__(self, edge_index, node_features, edge_features, node_labels, edge_labels):
        """
        Initialize a Graph instance.

        :param edge_index: A numpy array of shape (2, num_edges) representing the edges as pairs of node indices.
        :param node_features: A numpy array of node features with shape (num_nodes, num_node_features).
        :param edge_features: A numpy array of edge features with shape (num_edges, num_edge_features).
        :param node_labels: A numpy array of node labels with shape (num_nodes,).
        :param edge_labels: A numpy array of edge labels with shape (num_edges,).
        """
        self.node_features = np.array(node_features)
        self.edge_features = np.array(edge_features)
        self.edge_index = np.array(edge_index)
        self.node_labels = np.array(node_labels)
        self.edge_labels = np.array(edge_labels)

    def __repr__(self):
        return (f"Graph(node_features={self.node_features.shape}, "
                f"edge_features={self.edge_features.shape}, "
                f"edge_index={self.edge_index.shape}, "
                f"node_labels={self.node_labels.shape if self.node_labels is not None else None}, "
                f"edge_labels={self.edge_labels.shape if self.edge_labels is not None else None})")

# %%
def loadGraph(directory):
    """This function loads one specific node dataset and performs normalization
    Args:
        directory: path that contains the files truemult.node and truemult.edge
    
    """
    # ===================== NODES =====================
    raw_node = np.loadtxt(directory + "/truemult.node", dtype=str, delimiter='\t')

    # Extract the columns
    node_id  = raw_node[:,0].astype(int)
    node_cov = raw_node[:,1].astype(np.float32)
    node_len = raw_node[:,2].astype(np.float32)
    node_cg  = raw_node[:,3].astype(np.float32)
    node_tm  = raw_node[:,4].astype(int)

    # Compute the median coverage for nodes with true multiplicity = 1
    median_node_cov = np.median(node_cov[node_tm == 1])

    # Preprocess node data for training
    node_tm[node_tm > max_mult] = max_mult   # max_mult + 1 classes
    node_cov = node_cov / median_node_cov    # normalize
    node_len = np.log(node_len)              # log-transform node length

    # reorder and expand features such that they appear in node order
    # [-1, 1, -2, 2, -3, 3, ... -N, N] with N = #nodes
    node_cov_ordered = np.empty_like(node_cov)
    node_cov_ordered[node_id - 1] = node_cov
    node_cov_ordered = node_cov_ordered.repeat(2)

    node_cg_ordered = np.empty_like(node_cg)
    node_cg_ordered[node_id - 1] = node_cg
    node_cg_ordered = node_cg_ordered.repeat(2)

    node_len_ordered = np.empty_like(node_len)
    node_len_ordered[node_id - 1] = node_len
    node_len_ordered = node_len_ordered.repeat(2)

    node_labels = np.empty_like(node_tm)
    node_labels[node_id - 1] = node_tm
    node_labels = node_labels.repeat(2)  

    # ===================== EDGES =====================  
    raw_edge = np.loadtxt(directory + "/truemult.edge", dtype=str, delimiter='\t')

    # Extract the columns
    edge_src = raw_edge[:,0].astype(int)
    edge_dst = raw_edge[:,1].astype(int)
    edge_cov = raw_edge[:,2].astype(np.float32)
    edge_cg  = raw_edge[:,3].astype(np.float32)
    edge_tm  = raw_edge[:,4].astype(int)

    # Compute the median coverage for edges with true multiplicity = 1
    median_edge_cov = np.median(edge_cov[edge_tm == 1]) # TODO what with unseen data without GT?
    
    # Preprocess edge data for training
    edge_tm[edge_tm > max_mult] = max_mult
    edge_cov = edge_cov / median_edge_cov

    # Expand edge features such that for every edge (src, dst), 
    # the reverse-complementary edge (-dst, -src) is added.
    # The indexes in edge_index point to entries in node_X_ordered
    edge_index = np.zeros(shape=(2, 2 * edge_src.size), dtype=int)
    for i in range(0, edge_src.size):
        if (edge_src[i] > 0):
            edge_index[0, 2*i]   =  2*edge_src[i]-1  # fwd edge
            edge_index[1, 2*i+1] =  2*edge_src[i]-2  # rc edge
        else:
            edge_index[0, 2*i]   = -2*edge_src[i]-2
            edge_index[1, 2*i+1] = -2*edge_src[i]-1        
        if (edge_dst[i] > 0):
            edge_index[1, 2*i]   =  2*edge_dst[i]-1  # fwd edge
            edge_index[0, 2*i+1] =  2*edge_dst[i]-2  # rc edge
        else:
            edge_index[1, 2*i]   = -2*edge_dst[i]-2
            edge_index[0, 2*i+1] = -2*edge_dst[i]-1

    edge_cov_ordered = edge_cov.repeat(2)
    edge_cg_ordered  = edge_cg.repeat(2)
    edge_labels      = edge_tm.repeat(2)

    node_features = np.column_stack((node_cov_ordered, node_cg_ordered, node_len_ordered))
    edge_features = np.column_stack((edge_cov_ordered, edge_cg_ordered))

    print(f"Loaded graph with {node_id.size} nodes and {edge_src.size} edges") 
    #print(f"Median coverage is {median_node_cov:.2f} for nodes and {median_edge_cov:.2f} for edges")
    
    return Graph(edge_index, node_features, edge_features, node_labels, edge_labels)

# %%
def visualizeCoverage(features, labels):
    """This function visualizes the coverage distribution of nodes or edges.

    Args:
        features: either node_features or edge_features as obtained by loadGraph
        labels: either node_labels or edge_labels as obtained by loadGraph
    """

    # Undo the .repeat(2) in labels and features - visualize only 
    features = features[::2]
    labels = labels[::2]

    cov0 = features[labels == 0, 0]
    cov1 = features[labels == 1, 0]
    cov2 = features[labels == 2, 0]
    cov3 = features[labels == 3, 0]
    cov4 = features[labels == 4, 0]
    cov5 = features[labels == 5, 0]

    print(f"Number of data points with coverage 0 : {cov0.size}")
    print(f"Number of data points with coverage 1 : {cov1.size}")
    print(f"Number of data points with coverage 2 : {cov2.size}")
    print(f"Number of data points with coverage 3 : {cov3.size}")
    print(f"Number of data points with coverage 4 : {cov4.size}")
    print(f"Number of data points with coverage 5+: {cov5.size}")

    boxData = [cov0, cov1, cov2, cov3, cov4, cov5]

    fig, ax = plt.subplots()

    # Create a boxplot
    ax.boxplot(boxData)
    plt.xticks([1, 2, 3, 4, 5, 6], ['cov 0', 'cov 1', 'cov 2', 'cov 3', 'cov4', 'cov5'])
    plt.ylim(0, 10)
    yticks = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    dummy = plt.yticks(yticks)

    # Add horizontal lines at each y-tick
    for y in yticks:
        ax.axhline(y, color='gray', linestyle='dashed', linewidth=0.5)

# %%
def get_nodes_within_depth(G, start_node, depth):
    # Use BFS to find all nodes within the specified depth, considering both directions
    visited = set()
    queue = deque([(start_node, 0)])

    while queue:
        node, current_depth = queue.popleft()
        if current_depth > depth:
            continue
        if node not in visited:
            visited.add(node)
            # Include both successors and predecessors
            for neighbor in G.successors(node):
                if neighbor not in visited:
                    queue.append((neighbor, current_depth + 1))
            for neighbor in G.predecessors(node):
                if neighbor not in visited:
                    queue.append((neighbor, current_depth + 1))
    return visited

def visualizeSubgraph(start_node, depth, graph):
    """This function visualizes a subgraph

    Args:
        start_node: the central node of the subgraph (in Detox naming scheme)
        depth: the maximum path length between the central_node and any visualized node
        graph: as obtained by loadGraph
    """
    # Create a directed NetworkX graph from the edge index
    G = nx.DiGraph()
    G.add_edges_from(graph.edge_index.T)

    # Increase figure size
    plt.figure(figsize=(12, 4))

    # Convert Detox naming scheme to index in node_features / node_labels
    start_node = 2*start_node-1 if start_node > 0 else -2*start_node-2

    # Get nodes within the specified depth
    nodes_within_depth = get_nodes_within_depth(G, start_node, depth)

    # Build the node labels dictionary
    node_text = {}
    for node in nodes_within_depth:
        node_text[node] = str(-(node // 2 + 1) if node % 2 == 0 else (node + 1) // 2) #+ " (" + str(node_labels[node]) + ")"

    nx.set_node_attributes(G, node_text, 'label')
   # nx.set_edge_attributes(G, edge_labels, 'label')

    # Extract the subgraph
    subgraph = G.subgraph(nodes_within_depth)

    # Visualize the subgraph with arrows on the edges
    pos = nx.kamada_kawai_layout(subgraph)
    pos = nx.spring_layout(subgraph, pos=pos, iterations=3)
    nx.draw(subgraph, pos, with_labels=True, labels=nx.get_node_attributes(subgraph, 'label'),
            node_color='lightblue', edge_color='gray', node_size=1000, arrows=True)

    # Draw edge labels
    #nx.draw_networkx_edge_labels(subgraph, pos, edge_labels=nx.get_edge_attributes(subgraph, 'label'))

    plt.show()

# %%
class TGraph:
    def __init__(self, graph, device, train_frac = 0.9):
        """
        TGraph is identical to Graph except that it uses tensors and computes tensors
        to split the graph in a training and validation set
        This constructor accepts a Graph object and train_frac, the percentage
        of the number of nodes and edges that should be used for training
        """
        # fixed seed for PyTorch
        torch.manual_seed(seed)
        self.node_features = torch.tensor(graph.node_features, dtype=torch.float, device=device)
        self.edge_features = torch.tensor(graph.edge_features, dtype=torch.float, device=device)
        self.edge_index = torch.tensor(graph.edge_index, dtype=torch.long, device=device)
        self.node_labels = torch.tensor(graph.node_labels, dtype=torch.long, device=device)
        self.edge_labels = torch.tensor(graph.edge_labels, dtype=torch.long, device=device)

        # split the data in a training and validation set
        node_train_size = int(len(self.node_labels) * train_frac)                                  
        node_indices = torch.randperm(len(self.node_labels))
        self.node_train_indices = node_indices[:node_train_size]
        self.node_val_indices = node_indices[node_train_size:]
        self.node_train_labels = self.node_labels[self.node_train_indices]
        self.node_val_labels = self.node_labels[self.node_val_indices]

        edge_train_size = int(len(self.edge_labels) * train_frac)                                  
        edge_indices = torch.randperm(len(self.edge_labels))
        self.edge_train_indices = edge_indices[:edge_train_size]
        self.edge_val_indices = edge_indices[edge_train_size:]
        self.edge_train_labels = self.edge_labels[self.edge_train_indices]
        self.edge_val_labels = self.edge_labels[self.edge_val_indices]      

    def __repr__(self):
        return (f"TGraph(node_features={self.node_features.shape}, "
                f"edge_features={self.edge_features.shape}, "
                f"edge_index={self.edge_index.shape}, "
                f"node_labels={self.node_labels.shape}, "
                f"edge_labels={self.edge_labels.shape})")

class TGraph_kfold:
    def __init__(self, graph, device, seed, n_splits=10):
        """
        TGraph_kfold is identical to Graph except that it uses tensors and computes tensors
        to split the graph in a training and validation set
        This constructor accepts a Graph object and n_splits, the number of folds of splits 
        of nodes and edges that should be used for training
        """
        # fixed seed for PyTorch
        torch.manual_seed(seed)
        self.node_features = torch.tensor(graph.node_features, dtype=torch.float, device=device)
        self.edge_features = torch.tensor(graph.edge_features, dtype=torch.float, device=device)
        self.edge_index = torch.tensor(graph.edge_index, dtype=torch.long, device=device)
        self.node_labels = torch.tensor(graph.node_labels, dtype=torch.long, device=device)
        self.edge_labels = torch.tensor(graph.edge_labels, dtype=torch.long, device=device)

        node_even_indices = [2*i for i in range(int(len(graph.node_labels)/2))]
        edge_even_indices = [2*i for i in range(int(len(graph.edge_labels)/2))]
        
        # split the data in a training and validation set
        self.node_train_indices = {}
        self.node_val_indices = {}
        self.skf = StratifiedKFold(n_splits, shuffle=True, random_state=seed)
        for fold, (train_idx, val_idx) in enumerate(self.skf.split(self.node_features[node_even_indices].cpu(), self.node_labels[node_even_indices].cpu())):
            node_train_idx = [ti*2 for ti in train_idx] + [ti*2+1 for ti in train_idx]
            node_val_idx = [vi*2 for vi in val_idx] + [vi*2+1 for vi in val_idx]
            self.node_train_indices[fold] = node_train_idx
            self.node_val_indices[fold] = node_val_idx
        self.node_train_labels = {fold: self.node_labels[nti] for fold, nti in self.node_train_indices.items()}
        self.node_val_labels = {fold: self.node_labels[nvi] for fold, nvi in self.node_val_indices.items()}

        self.edge_train_indices = {}
        self.edge_val_indices = {}
        self.skf = StratifiedKFold(n_splits, shuffle=True, random_state=seed)
        for fold, (train_idx, val_idx) in enumerate(self.skf.split(self.edge_features[edge_even_indices].cpu(), self.edge_labels[edge_even_indices].cpu())):
            edge_train_idx = [ti*2 for ti in train_idx] + [ti*2+1 for ti in train_idx]
            edge_val_idx = [vi*2 for vi in val_idx] + [vi*2+1 for vi in val_idx]
            self.edge_train_indices[fold] = edge_train_idx
            self.edge_val_indices[fold] = edge_val_idx
        self.edge_train_labels = {fold: self.edge_labels[eti] for fold, eti in self.edge_train_indices.items()}
        self.edge_val_labels = {fold: self.edge_labels[evi] for fold, evi in self.edge_val_indices.items()}      

    def __repr__(self):
        return (f"TGraph_kfold(node_features={self.node_features.shape}, "
                f"edge_features={self.edge_features.shape}, "
                f"edge_index={self.edge_index.shape}, "
                f"node_labels={self.node_labels.shape}, "
                f"edge_labels={self.edge_labels.shape})")

# %%
class MLPnet(nn.Module):
    def __init__(self, num_node_features, num_edge_features, num_node_classes, num_edge_classes):
        super(MLPnet, self).__init__()
        self.nfc1 = nn.Linear(num_node_features, 16)
        self.nfc2 = nn.Linear(16, 8)
        self.nfc3 = nn.Linear(8, num_node_classes)
        self.efc1 = nn.Linear(num_edge_features, 16)
        self.efc2 = nn.Linear(16, 8)
        self.efc3 = nn.Linear(8, num_edge_classes)

    def forward(self, node_features, edge_index, edge_features):
        xn = F.relu(self.nfc1(node_features))
        xn = F.relu(self.nfc2(xn))
        xn = self.nfc3(xn)
        xe = F.relu(self.efc1(edge_features))
        xe = F.relu(self.efc2(xe))
        xe = self.efc3(xe)
        return xn, xe

# Simon's model adapted such that we can support everything with evidence from literature
class MessagePassingGNN_bis(torch.nn.Module):
    """
    Message-passing graph neural network for node and edge classification.

    Args:
        node_in_channels (int): Number of input features for each node.
        edge_in_channels (int): Number of input features for each edge.
        hidden_channels (int): Dimensionality of hidden representations.
        num_node_classes (int): Number of output classes for node classification.
        num_edge_classes (int): Number of output classes for edge classification.
        use_cg (bool): use CG-percentage as input feature
        use_len (bool): use number of k-mers in a node as input feature
        use_highway (str): which, if any, highway-gating to use. 
                            ('identity': no highway gates,
                            'sigmoid': highway gate with sigmoid,
                            'sigmoid_with_update': highway gate with sigmoid,
                                include both old and updated message in gating computation)
        feature_in_location (str): Where to use features other than coverage. 
                                    ('initial': as additional input,
                                    'update': as external feature in message updates,
                                    'end': at the very end when all messages have been passed)
        num_propagations (int): Number of message passing layers.
        activation (callable): Activation function to apply (e.g., torch.nn.functional.elu).
    """

    def __init__(self, node_in_channels, edge_in_channels, hidden_channels, 
                 num_node_classes, num_edge_classes, use_cg=True, use_len=False, use_highway="sigmoid_with_update",
                 feature_in_location = 'end', num_propagations=7, activation=F.elu):
        super(MessagePassingGNN_bis, self).__init__()

        if feature_in_location == 'initial':
            assert((use_cg + use_len + 1) == node_in_channels)
            assert((use_cg + 1) == edge_in_channels)
        self.feature_in = feature_in_location
        self.use_cg = use_cg
        self.use_len = use_len
        self.use_gate = ('sigmoid' in use_highway)
        self.gate_type = use_highway
        self.nfeat_node = 1 + use_cg + use_len
        self.nfeat_edge = 1 + use_cg

        self.num_propagations = num_propagations

        # TODO why a single hidden channels parameter that is the same in all linear layers used
        self.node_initial_mlp = Linear(node_in_channels, hidden_channels)
        self.edge_initial_mlp = Linear(edge_in_channels, hidden_channels)
        self.node_initial_ln = torch.nn.LayerNorm(hidden_channels)
        self.edge_initial_ln = torch.nn.LayerNorm(hidden_channels)
        self.message_mlp = Linear(hidden_channels * 2, hidden_channels)

        node_update_in_channels = hidden_channels + 4 * hidden_channels
        if self.feature_in == 'update':
            node_update_in_channels = node_update_in_channels + use_cg + use_len
        self.node_mlps = torch.nn.ModuleList([
            Linear(node_update_in_channels, hidden_channels) for _ in range(num_propagations)
        ])
        self.node_lns = torch.nn.ModuleList([
            torch.nn.LayerNorm(hidden_channels) for _ in range(num_propagations)
        ])

        edge_update_in_channels = 2 * hidden_channels + hidden_channels
        if self.feature_in == 'update':
            edge_update_in_channels = edge_update_in_channels + use_cg
        self.edge_mlps = torch.nn.ModuleList([
            Linear(edge_update_in_channels, hidden_channels) for _ in range(num_propagations)
        ])

        self.edge_lns = torch.nn.ModuleList([
            torch.nn.LayerNorm(hidden_channels) for _ in range(num_propagations)
        ])
        if self.use_gate:
            assert((self.gate_type in ["sigmoid", "sigmoid_with_update"]))
            if self.gate_type == 'sigmoid':
                print("using sigmoid gate on prev. h")
                self.node_gate_mlps = torch.nn.ModuleList([
                    Linear(hidden_channels, hidden_channels) for _ in range(num_propagations)
                ])
                self.edge_gate_mlps = torch.nn.ModuleList([
                    Linear(hidden_channels, hidden_channels) for _ in range(num_propagations)
                ])
            elif self.gate_type == 'sigmoid_with_update':
                print("using sigmoid gate on prev. h and update m")
                self.node_gate_mlps = torch.nn.ModuleList([
                    Linear(2*hidden_channels, hidden_channels) for _ in range(num_propagations)
                ])
                self.edge_gate_mlps = torch.nn.ModuleList([
                    Linear(2*hidden_channels, hidden_channels) for _ in range(num_propagations)
                ])
        else:
            self.node_gate_mlps = torch.nn.ModuleList([
                torch.nn.Identity() for _ in range(num_propagations)
            ])
            self.edge_gate_mlps = torch.nn.ModuleList([
                torch.nn.Identity() for _ in range(num_propagations)
            ])
            
        node_class_in_channels = hidden_channels
        edge_class_in_channels = hidden_channels
        if (self.feature_in == 'end') or (self.feature_in == 'end_mlp'):
            node_class_in_channels = node_class_in_channels + use_cg + use_len
            edge_class_in_channels = edge_class_in_channels + use_cg
        
        if (self.feature_in == 'end_mlp'):
            self.node_classifier = nn.Sequential(Linear(node_class_in_channels, 2*hidden_channels),
                                                 nn.ELU(),
                                                 Linear(2*hidden_channels, num_node_classes))
            self.edge_classifier = nn.Sequential(Linear(edge_class_in_channels, 2*hidden_channels),
                                                 nn.ELU(),
                                                 Linear(2*hidden_channels, num_edge_classes))
        else:
            self.node_classifier = Linear(node_class_in_channels, num_node_classes)
            self.edge_classifier = Linear(edge_class_in_channels, num_edge_classes)
        self.activation = activation

    def forward(self, x, edge_index, edge_attr):
        """
        Run forward propagation through the MPNN model.

        Args:
            x (Tensor): Node input features of shape (num_nodes, 3).
            edge_index (Tensor): Edge indices (2, num_edges).
            edge_attr (Tensor): Edge input features of shape (num_edges, 2).

        Returns:
            tuple[Tensor, Tensor]:
                - node_out: Node classification logits of shape (num_nodes, num_node_classes).
                - edge_out: Edge classification logits of shape (num_edges, num_edge_classes).
        """
        cov = x[:, 0].unsqueeze(1)
        cg  = x[:, 1].unsqueeze(1)
        length = x[:, 2].unsqueeze(1) ## Unused (see length test in ablation study Simon)
        
        node_in = cov
        if self.feature_in == 'initial':
            if self.use_cg:
                node_in = torch.cat([node_in, cg], dim=1)
            if self.use_len:
                node_in = torch.cat([node_in, length], dim=1)

        cov_edge = edge_attr[:, 0].unsqueeze(1)
        cg_edge = edge_attr[:, 1].unsqueeze(1)

        edge_in = cov_edge
        if self.feature_in == 'initial':
            if self.use_cg:
                edge_in = torch.cat([edge_in, cg_edge], dim=1)

        # "obtain a representative feature vector appropriate for further processing"
        node_features = self.node_initial_ln(self.activation(self.node_initial_mlp(node_in)))
        edge_features = self.edge_initial_ln(self.activation(self.edge_initial_mlp(edge_in)))

        for node_mlp, node_ln, edge_mlp, edge_ln, node_gate, edge_gate in zip(self.node_mlps, self.node_lns, self.edge_mlps, self.edge_lns, self.node_gate_mlps, self.edge_gate_mlps):
            # Get node messages from incoming and outgoing neighbours
            combined_messages = self.aggregate_messages(node_features, edge_features, edge_index)
            # compute node update using prev t features and messages
            update_node_in = torch.cat([node_features, combined_messages], dim=1)
            if self.feature_in == 'update':
                if self.use_cg:
                    update_node_in = torch.cat([update_node_in, cg], dim=1)
                if self.use_len:
                    update_node_in = torch.cat([update_node_in, length], dim=1)
            node_update = node_ln(self.activation(node_mlp(update_node_in)))
            # highway gate TODO: do we need a gate?
            ## based on Highway GCN updates https://arxiv.org/pdf/1804.08049, https://www.sciencedirect.com/science/article/pii/S2666651021000012 https://en.wikipedia.org/wiki/Highway_network
            ## TODO: test input of node_update into node_transform_gate
            if self.use_gate:
                if self.gate_type == 'sigmoid':
                    node_transform_gate = F.sigmoid(node_gate(node_features)) # resemble highway network formulation better
                elif self.gate_type == 'sigmoid_with_update':
                    node_transform_gate = F.sigmoid(node_gate(torch.cat([node_features, node_update], dim=1)))
                else:
                    node_transform_gate = F.sigmoid(node_gate(node_features))
                node_features = node_transform_gate  * node_update + (1 - node_transform_gate) * node_features
            else: 
                node_features = node_update
            # update edges
            source_features = node_features[edge_index[0]] # TODO: edge features use updated node features, is not used in Gilmer, but is used as such in https://arxiv.org/pdf/1806.03146
            target_features = node_features[edge_index[1]]
            update_edge_in = torch.cat([edge_features, source_features, target_features], dim=1)
            if self.feature_in == 'update':
                if self.use_cg:
                    update_edge_in = torch.cat([update_edge_in, cg_edge], dim=1)
            edge_update = edge_ln(self.activation(edge_mlp(update_edge_in)))
            if self.use_gate:
                if self.gate_type == 'sigmoid':
                    edge_transform_gate = F.sigmoid(edge_gate(edge_features))
                elif self.gate_type == 'sigmoid_with_update':
                    edge_transform_gate = F.sigmoid(edge_gate(torch.cat([edge_features, edge_update], dim=1)))
                else:
                    edge_transform_gate = F.sigmoid(edge_gate(edge_features))
                edge_features = edge_transform_gate *  edge_update + (1 - edge_transform_gate) * edge_features
            else:
                edge_features = edge_update

        # TODO add cg and length as feature in node/edge update or add into initial mlp as feature
        if (self.feature_in == 'end') or (self.feature_in == 'end_mlp'):
            if self.use_cg:
                node_features = torch.cat([node_features, cg], dim=1)
                edge_features = torch.cat([edge_features, cg_edge], dim=1)
            if self.use_len:
                node_features = torch.cat([node_features, length], dim=1)
        
        node_out = self.node_classifier(node_features)
        edge_out = self.edge_classifier(edge_features)

        return node_out, edge_out

    def aggregate_messages(self, node_features, edge_features, edge_index):
        """
        Aggregate messages from source and target nodes for each edge.

        Args:
            node_features (Tensor): Node hidden states of shape (num_nodes, hidden_dim).
            edge_features (Tensor): Edge hidden states of shape (num_edges, hidden_dim).
            edge_index (Tensor): Edge indices (2, num_edges).

        Returns:
            Tensor: Concatenated incoming and outgoing message aggregates per node,
                    shape (num_nodes, 4 * hidden_dim).
        """
        i = edge_index[0].to(torch.int64)
        j = edge_index[1].to(torch.int64)

        source_feat = node_features[i]
        target_feat = node_features[j]
                
        messages_from_source = torch.cat([edge_features, source_feat], dim=1)
        messages_from_target = torch.cat([edge_features, target_feat], dim=1)

        # print(j.unsqueeze(-1).expand(-1, messages_from_source.size(1)))

        num_nodes = node_features.size(0)

        # print(j.unsqueeze(-1).expand(-1, messages_from_source.size(1)).size(), messages_from_source.size())

        incoming = torch.zeros((num_nodes, messages_from_source.size(1)), device=messages_from_source.device)
        incoming.scatter_add_(0, j.unsqueeze(-1).expand(-1, messages_from_source.size(1)), messages_from_source)
        outgoing = torch.zeros((num_nodes, messages_from_target.size(1)), device=messages_from_target.device)
        outgoing.scatter_add_(0, i.unsqueeze(-1).expand(-1, messages_from_target.size(1)), messages_from_target)

        combined = torch.cat([incoming, outgoing], dim=1) # TODO message could also be governed by a neural network instead of a simple sum
        # print(combined.size())
        return combined

# %%
# TODO add tensorboard tracking
def trainModel(model, data, optimizer, node_criterion, edge_criterion, device='cpu', swriter=None, num_epochs = 300, patience = 10):
    best_val_loss = float('inf')
    best_model_weights = None
    epochs_without_improvement = 0

    for epoch in range(num_epochs):
        # Shuffle the graphs at the beginning of each epoch
        random.shuffle(data)
        
        # Training phase
        model.train()
        total_train_loss = 0
        for g in data:    # batch = a single graph
            optimizer.zero_grad()
            node_out, edge_out = model.forward(g.node_features, g.edge_index, g.edge_features)
            node_loss = node_criterion(node_out[g.node_train_indices], g.node_train_labels)
            edge_loss = edge_criterion(edge_out[g.edge_train_indices], g.edge_train_labels)
            train_batch_loss = 0.5 * node_loss + 0.5 * edge_loss
            train_batch_loss.backward()
            optimizer.step()
            total_train_loss += train_batch_loss.item()
        avg_train_loss = total_train_loss / len(data)
        swriter.add_scalar('loss/training', avg_train_loss, epoch)

        # Validation phase
        model.eval()
        total_val_loss = 0
        with torch.no_grad():
            for g in data:
                node_out, edge_out = model.forward(g.node_features, g.edge_index, g.edge_features)
                node_loss = node_criterion(node_out[g.node_val_indices], g.node_val_labels)
                edge_loss = edge_criterion(edge_out[g.edge_val_indices], g.edge_val_labels)
                val_batch_loss = 0.5 * node_loss + 0.5 * edge_loss
                total_val_loss += val_batch_loss.item()
        avg_val_loss = total_val_loss / len(data)
        swriter.add_scalar('loss/validation', avg_val_loss, epoch)

        # Print training and validation loss
        if epoch % 10 == 0:
            print(f"Epoch: {epoch}; training loss: {avg_train_loss}; validation Loss: {avg_val_loss}")
        
        # Save the best model weights based on validation loss
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            best_model_weights = copy.deepcopy(model.state_dict())
            epochs_without_improvement = 0
        else:
            epochs_without_improvement += 1


        # Early stopping check
        if epochs_without_improvement >= patience:
            print(f"Early stopping criterion triggered after {epoch} epochs.")
            break

    if best_model_weights is not None:
        model.load_state_dict(best_model_weights)
    print(f"Best validation loss: {best_val_loss}")
    metric_dict = {'best validation loss': best_val_loss, 'number of epochs trained': epoch,
                   'final training loss': avg_train_loss, 'final validation loss': avg_val_loss}
    return metric_dict

def trainModel_kfold(model, data, optimizer, node_criterion, edge_criterion, fold_idx,
               device='cpu', swriter=None, num_epochs = 300, patience = 10):
    best_val_loss = float('inf')
    best_model_weights = None
    epochs_without_improvement = 0

    for epoch in range(num_epochs):
        # Shuffle the graphs at the beginning of each epoch
        random.shuffle(data)
        
        # Training phase
        model.train()
        total_train_loss = 0
        for g in data:    # batch = a single graph
            optimizer.zero_grad()
            node_out, edge_out = model.forward(g.node_features, g.edge_index, g.edge_features)
            node_loss = node_criterion(node_out[g.node_train_indices[fold_idx]], g.node_train_labels[fold_idx])
            edge_loss = edge_criterion(edge_out[g.edge_train_indices[fold_idx]], g.edge_train_labels[fold_idx])
            train_batch_loss = 0.5 * node_loss + 0.5 * edge_loss
            train_batch_loss.backward()
            optimizer.step()
            total_train_loss += train_batch_loss.item()
        avg_train_loss = total_train_loss / len(data)
        swriter.add_scalar('loss/training', avg_train_loss, epoch)

        # Validation phase
        model.eval()
        total_val_loss = 0
        with torch.no_grad():
            for g in data:
                node_out, edge_out = model.forward(g.node_features, g.edge_index, g.edge_features)
                node_loss = node_criterion(node_out[g.node_val_indices[fold_idx]], g.node_val_labels[fold_idx])
                edge_loss = edge_criterion(edge_out[g.edge_val_indices[fold_idx]], g.edge_val_labels[fold_idx])
                val_batch_loss = 0.5 * node_loss + 0.5 * edge_loss
                total_val_loss += val_batch_loss.item()
        avg_val_loss = total_val_loss / len(data)
        swriter.add_scalar('loss/validation', avg_val_loss, epoch)

        # Print training and validation loss
        if epoch % 10 == 0:
            print(f"Epoch: {epoch}; training loss: {avg_train_loss}; validation Loss: {avg_val_loss}")
        
        # Save the best model weights based on validation loss
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            best_model_weights = copy.deepcopy(model.state_dict())
            epochs_without_improvement = 0
        else:
            epochs_without_improvement += 1


        # Early stopping check
        if epochs_without_improvement >= patience:
            print(f"Early stopping criterion triggered after {epoch} epochs.")
            break

    if best_model_weights is not None:
        model.load_state_dict(best_model_weights)
    print(f"Best validation loss: {best_val_loss}")
    metric_dict = {'best validation loss': best_val_loss, 'number of epochs trained': epoch,
                   'final training loss': avg_train_loss, 'final validation loss': avg_val_loss}
    return metric_dict

# %%
def evaluateModelAccuracy(model, test_data, node_criterion, edge_criterion, device='cpu'):
    total_nodes = 0
    correct_nodes = 0
    total_edges = 0
    correct_edges = 0

    model.eval()
    total_test_loss = 0
    with torch.no_grad():
        for g in test_data:
            node_out, edge_out = model.forward(g.node_features, g.edge_index, g.edge_features)
            node_loss = node_criterion(node_out, g.node_labels)
            edge_loss = edge_criterion(edge_out, g.edge_labels)
            test_batch_loss = 0.5 * node_loss + 0.5 * edge_loss
            total_test_loss += test_batch_loss.item()
    
            # count correct nodes
            _, predicted_nodes = torch.max(node_out.data, 1)
            total_nodes += g.node_labels.size(0)
            correct_nodes += (predicted_nodes == g.node_labels).sum().item()

            # count correct edges
            _, predicted_edges = torch.max(edge_out.data, 1)
            total_edges += g.edge_labels.size(0)
            correct_edges += (predicted_edges == g.edge_labels).sum().item()

    avg_test_loss = total_test_loss / len(test_data)

    print(f"Average test loss: {avg_test_loss:.4f}")
    print(f'Node Test Accuracy: {correct_nodes / total_nodes * 100:.2f}%')
    print(f'Edge Test Accuracy: {correct_edges / total_edges * 100:.2f}%')
    evaluation_dict = {'average test loss': avg_test_loss, 'node test accuracy': correct_nodes / total_nodes * 100, 'edge test accuracy': correct_edges / total_edges * 100}
    return evaluation_dict

# %%
def makeConfusionMatrix(model, test_data, save_loc='.', save_suffix=""):
    node_confusion_matrix = np.zeros((max_mult + 1, max_mult + 1))
    edge_confusion_matrix = np.zeros((max_mult + 1, max_mult + 1))

    model.eval()
    with torch.no_grad():
        for g in test_data: 
            node_out, edge_out = model.forward(g.node_features, g.edge_index, g.edge_features)
    
            # count correct nodes
            _, predicted_nodes = torch.max(node_out.data, 1)
            for l, p in zip(g.node_labels.view(-1), predicted_nodes.view(-1)):
                node_confusion_matrix[l, p] += 1

            # count correct edges
            _, predicted_edges = torch.max(edge_out.data, 1)
            for l, p in zip(g.edge_labels.view(-1), predicted_edges.view(-1)):
                edge_confusion_matrix[l, p] += 1

    # Normalize confusion matrix by row (true label percentages)
    node_conf_matrix_normalized = node_confusion_matrix / node_confusion_matrix.sum(axis=1, keepdims=True)
    edge_conf_matrix_normalized = edge_confusion_matrix / edge_confusion_matrix.sum(axis=1, keepdims=True)
    ytick_labels = list(range(node_confusion_matrix.shape[0]))

    plt.figure(figsize=(8, 6))
    sns.heatmap(
        node_conf_matrix_normalized,
        annot=True,  
        fmt=".2f", 
        cmap="Greys",  
        xticklabels=range(node_confusion_matrix.shape[1]),
        yticklabels=ytick_labels,
        cbar_kws={'label': 'Percentage'},
        square=True, 
        vmin=0, vmax=1 
    )
    plt.xlabel("Predicted Classes")
    plt.ylabel("True Classes")
    plt.title("Normalized Confusion Matrix Heatmap, Nodes")
    plt.savefig(f"{save_loc}/nodes_confusion{save_suffix}.png")
    # plt.show()
    plt.close()

    plt.figure(figsize=(8, 6))
    sns.heatmap(
        edge_conf_matrix_normalized,
        annot=True,
        fmt=".2f",
        cmap="Greys",
        xticklabels=range(edge_confusion_matrix.shape[1]),
        yticklabels=ytick_labels,
        cbar_kws={'label': 'Percentage'},
        square=True,
        vmin=0, vmax=1
    )
    plt.xlabel("Predicted Classes")
    plt.ylabel("True Classes")
    plt.title("Normalized Confusion Matrix Heatmap, Edge")
    plt.savefig(f"{save_loc}/edges_confusion{save_suffix}.png")
    # plt.show()
    plt.close()

    print(node_confusion_matrix)
    print(edge_confusion_matrix)

# %%
if __name__== "__main__":
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print('Using device:', device)

    parser = argparse.ArgumentParser()

    parser.add_argument('-s', '--seed', type=int, help='seed', default = 42)
    parser.add_argument('-f', '--featureloc', type=str, help='where to input features other than coverage ("initial", "update", "end", "end_mlp")', default = 'end')
    parser.add_argument('-g', '--gate', type=str, help='gating strategy ("identity" (=no gate), "sigmoid" (=highway gate), "sigmoid_with_update")', default = 'sigmoid_with_update')
    parser.add_argument('-d', '--nlayers', type=int, help='number of mpnn layers to use', default = 7)
    parser.add_argument('-w', '--nhidden', type=int, help='number of nodes in the hidden layers of used MLPs', default = 8)
    parser.add_argument('--nocg', default=False, action=argparse.BooleanOptionalAction,
                        help="Use cg percentage as extra feature")
    parser.add_argument('--uselen', default=False, action=argparse.BooleanOptionalAction,
                        help="Use node length as extra feature")


    args = parser.parse_args()
    seed = args.seed
    input_addit_features_at = args.featureloc # 'initial', 'update', 'end', 'end_mlp'
    use_gate = args.gate # 'identity' (=no gate), 'sigmoid' (=highway gate), 'sigmoid_with_update'
    use_cg = (not args.nocg)
    use_len = args.uselen
    num_layers = args.nlayers
    num_hidden = args.nhidden

    if input_addit_features_at == 'initial':
        node_in_c = 1 + use_cg + use_len
        edge_in_c = 1 + use_cg
    else:
        node_in_c = 1
        edge_in_c = 1

    max_mult = 5
    nfolds = 10
    
    numLayers = num_layers  # number of layers in the GNN
    tensorboard_dir_suffix = f'{datetime.datetime.now().strftime("%d%b%y_%H-%M")}_seed{seed}'
    print("Tracking training statistics in tensorboard dir.:",  tensorboard_dir_suffix)
    models_dir = f"models_{tensorboard_dir_suffix}"
    os.mkdir(models_dir)
    for f in range(nfolds):
        os.mkdir(f"{models_dir}/fold{f}")
    ### Hyperparameters ###
    max_epoch = 1000
    lr = 0.001
    label_smooth_param = 0.1
    n_hidden_chan = num_hidden
    

    hyperparams_mlp = {'max_epoch': max_epoch, 'learning_rate': lr, "label_smoothing": label_smooth_param,
                    "use_cg": True, "use_len": True, "max_multiplicity": max_mult}
    hyperparams_gnn = {'max_epoch': max_epoch, 'learning_rate': lr, "label_smoothing": label_smooth_param,
                    "use_cg": use_cg, "use_len": use_len, "max_multiplicity": max_mult,
                    'activation': 'ELU', 'n_mp_layers': numLayers, "n_hidden_features": n_hidden_chan,
                    'additional_features_input_location': input_addit_features_at, 'used_gate': use_gate}

    print(hyperparams_gnn)

    t_graphs = {}

    for dir in directories:
        graph = loadGraph(data_folder + dir + '/REAL/30/') ## TODO hard coded the location
        t_graphs[dir] = TGraph_kfold(graph, device, seed, n_splits=nfolds)
        
    for test_org in groups:
        bac_train_val = [gp for gp in groups if not gp in [test_org]]
        bac_test = [test_org]
    
        t_train_val_data = []
        t_test_data = []
        
        for dir in directories:
            if any(bac in dir for bac in bac_train_val):
                t_train_val_data.append(t_graphs[dir])
            if any(bac in dir for bac in bac_test):
                t_test_data.append(t_graphs[dir])
                
        print(f'Training/validating on {bac_train_val}, testing on {bac_test}')
    
        print(f'There are {len(t_train_val_data)} graphs in the training/validation dataset')
        print(f'There are {len(t_test_data)} graphs in the test dataset', flush=True)
        for f in range(nfolds):
            print(f"Training on fold {f}")
            swriter = SummaryWriter(log_dir= f"runs/mpnn/new_mpnn_kfold/{test_org}/{tensorboard_dir_suffix}/fold{f}")
            hyperparams_gnn.update({'testset': test_org, 'fold':f})
        
            node_criterion = nn.CrossEntropyLoss(label_smoothing=label_smooth_param)
            edge_criterion = nn.CrossEntropyLoss(label_smoothing=label_smooth_param)
        
            # train the GNN model
        
            # Set the seed for the random module (shuffling of graphs at the beginning of each epoch)
            random.seed(seed)
            # create a random number generator with fixed seed for reproducibility
            np.random.seed(seed)
            # fixed seed for PyTorch (random initialization of model weights)
            torch.manual_seed(seed)
            modelGNN = MessagePassingGNN_bis(node_in_channels=node_in_c, edge_in_channels=edge_in_c, hidden_channels=n_hidden_chan, 
                                             num_node_classes=max_mult + 1, num_edge_classes=max_mult + 1,
                                             use_cg= use_cg, use_len= use_len, use_highway=use_gate,
                                             feature_in_location=input_addit_features_at,
                                             num_propagations=numLayers, activation=torch.nn.functional.elu).to(device)
            # print(summary(modelGNN))
            optimizer = optim.Adam(modelGNN.parameters(), lr=lr,
                                   betas=(0.9, 0.999), eps=1e-08, weight_decay=0) # lr=0.01 in train_within_dataset
            mpnn_metrics = trainModel_kfold(modelGNN, t_train_val_data, optimizer, node_criterion, edge_criterion, f, swriter=swriter, num_epochs = max_epoch, patience = 20)
        
            # save the model
            scripted_model = torch.jit.script(modelGNN.to('cpu'))
            scripted_model.save(f'{models_dir}/fold{f}/model_GNN_{test_org}.pt')
        
            print("GNN model:")
            mpnn_eval = evaluateModelAccuracy(modelGNN.to(device), t_test_data, node_criterion, edge_criterion)
        
            mpnn_metrics.update(mpnn_eval)
            swriter.add_hparams(hyperparams_gnn, mpnn_metrics)
    
            makeConfusionMatrix(modelGNN, t_test_data, models_dir, f"_{tensorboard_dir_suffix}_{test_org}_{f}_GNN")




