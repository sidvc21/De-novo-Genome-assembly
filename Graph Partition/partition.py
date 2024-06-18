import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csgraph
from scipy.sparse.linalg import eigsh
from scipy.sparse.linalg import eigs


from scipy.sparse.linalg import eigsh, ArpackNoConvergence

#SPECTRAL PARTITIONING

def spectral_partitioning(G):
    # Construct the adjacency matrix
    A = nx.adjacency_matrix(G).astype(float)
    # Construct the degree matrix
    D = np.diag(np.array(A.sum(axis=1)).flatten())
    # Construct the Laplacian matrix
    L = D - A

    # Compute the eigenvalues and eigenvectors of the Laplacian matrix
    eigenvalues, eigenvectors = eigsh(L, k=2, which='SM')
    second_eigenvector = eigenvectors[:, 1]

    # Partition the nodes based on the second smallest eigenvector
    partition1 = [node for i, node in enumerate(G.nodes()) if second_eigenvector[i] >= 0]
    partition2 = [node for i, node in enumerate(G.nodes()) if second_eigenvector[i] < 0]

    return partition1, partition2

#RECURSIVE PARTITIONING
def partition_graph_recursive(G, max_nodes=30, current_subgraph_name="Subgraph", partition_level=0):
    if G.number_of_nodes() <= max_nodes:
        return [(current_subgraph_name, G)]

    # Perform spectral partitioning
    partition1, partition2 = spectral_partitioning(G)

    # Increment partition level
    partition_level += 1

    # Create subgraphs
    subgraph1 = G.subgraph(partition1).copy()
    subgraph2 = G.subgraph(partition2).copy()

    # Construct names for subgraphs
    subgraph_name_1 = f"{current_subgraph_name} Partition-{partition_level} 1"
    subgraph_name_2 = f"{current_subgraph_name} Partition-{partition_level} 2"

    # Recursively partition subgraphs
    subgraphs = []
    subgraphs.extend(partition_graph_recursive(subgraph1, max_nodes, subgraph_name_1, partition_level))
    subgraphs.extend(partition_graph_recursive(subgraph2, max_nodes, subgraph_name_2, partition_level))

    return subgraphs



#PRINTING SUBGRAPHS
def print_subgraphs(subgraphs):
    for idx, subgraph in enumerate(subgraphs):
        print(f"Subgraph {idx + 1}")
        print("Node 1\tNode 2\tEdge Weight")
        for u, v, data in subgraph.edges(data=True):
            print(f"{u}\t{v}\t{data['weight']}")
        print("\n")

#MAIN EXECUTION

import networkx as nx
import matplotlib.pyplot as plt

def visualize_de_bruijn(G):
    plt.figure(figsize=(8, 6))
    pos = nx.spring_layout(G)
    edge_labels = {(u, v): d['weight'] for u, v, d in G.edges(data=True)}
    nx.draw(G, pos, with_labels=True, node_size=1200, node_color='lightblue', font_size=12, font_weight='bold', edge_color='gray', arrowsize=20)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    plt.title('De Bruijn Graph Visualization')
    plt.show()

def read_fasta_file(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence)
    return sequences

#PRINTING SUBGRAPHS
def print_subgraphs(subgraphs):
    for idx, (subgraph_name, subgraph) in enumerate(subgraphs):
        num_nodes = subgraph.number_of_nodes()
        print(f"Subgraph {idx + 1} - Total Nodes: {num_nodes}")
        print("Node 1\tNode 2\tEdge Weight")
        for u, v, data in subgraph.edges(data=True):
            print(f"{u}\t{v}\t{data['weight']}")
        print("\n")

        # Visualize subgraph
        plt.figure(figsize=(8, 6))
        pos = nx.spring_layout(subgraph)
        edge_labels = {(u, v): d['weight'] for u, v, d in subgraph.edges(data=True)}
        nx.draw(subgraph, pos, with_labels=True, node_size=1200, node_color='lightblue', font_size=12, font_weight='bold', edge_color='gray', arrowsize=20)
        nx.draw_networkx_edge_labels(subgraph, pos, edge_labels=edge_labels)
        plt.title(f'Subgraph {idx + 1} Visualization')
        plt.show()
# Main Execution

file_path = r"/content/sars_cov2_100.fasta"
sequences = read_fasta_file(file_path)

k = 4  # Adjust the value of k as needed

G = de_bruijn(sequences, k)

visualize_de_bruijn(G)

original_name = "Original graph"
subgraphs = partition_graph_recursive(G, max_nodes=30, current_subgraph_name=original_name)

num_nodes = G.number_of_nodes()
print("Number of nodes in the graph:", num_nodes)

num_edges = G.number_of_edges()
print("Number of edges in the graph:", num_edges)

if num_nodes > 30:
    subgraphs = partition_graph_recursive(G, max_nodes=30)
    print_subgraphs(subgraphs)
else:
    print("No partitioning needed. The graph has 30 or fewer nodes.")
