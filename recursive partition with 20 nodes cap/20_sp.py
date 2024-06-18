import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
from Bio import SeqIO
import numpy as np
from sklearn.cluster import KMeans

# Function to generate k-mers from a sequence
def generate_kmers(sequence, k):
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    return kmers

# Function to construct a directed de Bruijn graph from a list of k-mers
def construct_de_bruijn_graph(kmers):
    graph = defaultdict(list)
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph[prefix].append(suffix)
    return graph

# Function to visualize the directed de Bruijn graph
def visualize_de_bruijn_graph(graph):
    G = nx.DiGraph()
    for node, neighbors in graph.items():
        for neighbor in neighbors:
            G.add_edge(node, neighbor)
    nx.draw(G, with_labels=True, node_size=700, node_color="skyblue", font_size=10, font_weight="bold")
    plt.title("Directed de Bruijn Graph")
    plt.show()

# Function to identify the number of incoming and outgoing edges for each node
def calculate_node_degrees(graph):
    in_degrees = {}
    out_degrees = {}
    for node, neighbors in graph.items():
        out_degrees[node] = len(neighbors)
        in_degrees[node] = sum(1 for kmer in all_kmers if kmer.endswith(node))
    return in_degrees, out_degrees

# Function to identify the node with the least number of incoming edges
def find_node_least_incoming_edges(graph):
    min_in_edges = float('inf')
    node_with_min_in_edges = None
    for node in graph.keys():
        in_edges = sum(1 for kmer in all_kmers if kmer.endswith(node))
        if in_edges < min_in_edges:
            min_in_edges = in_edges
            node_with_min_in_edges = node
    return node_with_min_in_edges, min_in_edges

# Function to identify the node with the least number of outgoing edges
def find_node_least_outgoing_edges(graph):
    min_out_edges = float('inf')
    node_with_min_out_edges = None
    for node, neighbors in graph.items():
        out_edges = len(neighbors)
        if out_edges < min_out_edges:
            min_out_edges = out_edges
            node_with_min_out_edges = node
    return node_with_min_out_edges, min_out_edges

# Read the reads from the FASTA file
reads = []
with open("/content/reads_random.fasta", "r") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        reads.append(str(record.seq))

# Generate k-mers from the reads
k = 4
all_kmers = []
for read in reads:
    all_kmers.extend(generate_kmers(read, k))

# Construct the directed de Bruijn graph
de_bruijn_graph = construct_de_bruijn_graph(all_kmers)

# Calculate the number of incoming and outgoing edges for each node
node_in_degrees, node_out_degrees = calculate_node_degrees(de_bruijn_graph)

# Find the node with the least number of incoming edges
node_least_incoming, min_incoming_edges = find_node_least_incoming_edges(de_bruijn_graph)
# Find the node with the least number of outgoing edges
node_least_outgoing, min_outgoing_edges = find_node_least_outgoing_edges(de_bruijn_graph)

if node_least_incoming == node_least_outgoing:
    graph_temp = de_bruijn_graph.copy()
    graph_temp.pop(node_least_outgoing)
    node_least_outgoing, min_outgoing_edges = find_node_least_outgoing_edges(graph_temp)

print("Node with least incoming edges:", node_least_incoming, "(", min_incoming_edges, "incoming edges)")
print("Node with least outgoing edges:", node_least_outgoing, "(", min_outgoing_edges, "outgoing edges)")

def count_edge_connections(graph):
    connections = []
    for node1 in graph:
        for node2 in graph:
            if node1 != node2:  # Exclude self-loops
                num_connections = sum(1 for neighbor in graph[node1] if neighbor == node2)
                if num_connections > 0:
                    connections.append((node1, node2, num_connections))
    return connections

edge_connections = count_edge_connections(de_bruijn_graph)

def visualize_de_bruijn_graph_with_weights(graph, edge_weights, title):
    G = nx.DiGraph()
    pos = nx.spring_layout(graph)
    nx.draw(graph, pos, with_labels=True, node_size=700, node_color="skyblue", font_size=10, font_weight="bold")
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_weights, font_color='red')
    plt.title(title)
    plt.show()

G = nx.DiGraph()
for connection in edge_connections:
    G.add_edge(connection[0], connection[1], weight=connection[2])

visualize_de_bruijn_graph_with_weights(G, nx.get_edge_attributes(G, 'weight'), "Directed de Bruijn Graph with Edge Weights")

# Function to convert de Bruijn graph into weighted adjacency matrix
def de_bruijn_to_weighted_adjacency(de_bruijn_graph):
    nodes = sorted(de_bruijn_graph.keys())
    n = len(nodes)
    adjacency_matrix = np.zeros((n, n))
    for i, node in enumerate(nodes):
        for neighbor in de_bruijn_graph[node]:
            j = nodes.index(neighbor)
            weight = de_bruijn_graph[node].count(neighbor)  # Weight based on overlap frequency
            adjacency_matrix[i, j] = weight
    return adjacency_matrix

# Convert de Bruijn graph to weighted adjacency matrix
weighted_adjacency_matrix = de_bruijn_to_weighted_adjacency(de_bruijn_graph)

# Print the weighted adjacency matrix
print("Weighted Adjacency Matrix:")
print(weighted_adjacency_matrix)

# Calculate the Laplacian matrix
n = weighted_adjacency_matrix.shape[0]
degree_matrix = np.diag(np.sum(weighted_adjacency_matrix, axis=1))
laplacian_matrix = degree_matrix - weighted_adjacency_matrix

# Calculate the eigenvectors corresponding to the smallest eigenvalues
eigenvalues, eigenvectors = np.linalg.eigh(laplacian_matrix)

# Find the index of the second smallest eigenvalue
second_smallest_eigenvalue_index = np.argsort(eigenvalues)[1]

# Check if the eigenvector corresponding to the second smallest eigenvalue contains both positive and negative values
eigenvector = eigenvectors[:, second_smallest_eigenvalue_index]
if np.any(eigenvector > 0) and np.any(eigenvector < 0):
    print("Eigenvector corresponding to the second smallest eigenvalue contains both positive and negative values.")
else:
    print("Eigenvector corresponding to the second smallest eigenvalue does not contain both positive and negative values.")

# Use KMeans clustering on the eigenvectors to partition the graph into two clusters
kmeans = KMeans(n_clusters=2, random_state=0).fit(eigenvectors[:, 1].reshape(-1, 1))
labels = kmeans.labels_

# Split nodes based on the clustering labels
nodes_partition1 = [node for node, label in zip(de_bruijn_graph.keys(), labels) if label == 0]
nodes_partition2 = [node for node, label in zip(de_bruijn_graph.keys(), labels) if label == 1]

# Count the total number of nodes in each partition
total_nodes_partition1 = len(nodes_partition1)
total_nodes_partition2 = len(nodes_partition2)

print("\nPartition 1 nodes:", nodes_partition1)
print("Total nodes in Partition 1:", total_nodes_partition1)
print("\nPartition 2 nodes:", nodes_partition2)
print("Total nodes in Partition 2:", total_nodes_partition2)
print(eigenvalues)
print(eigenvectors)

# Calculate the initial number of edges in the graph before partitioning
initial_edges_count = np.count_nonzero(weighted_adjacency_matrix)

# Print the initial number of edges in the graph
print("\nInitial number of edges in the graph before partitioning:", initial_edges_count)

# Calculate the number of edges cut off due to partitioning
edges_cut_off = 0
cut_off_edges_count = 0  # Counter for edges with non-zero weights
cut_off_edges_with_weight_gt_zero = 0  # Counter for edges with non-zero weights greater than zero
for node_index, node in enumerate(de_bruijn_graph.keys()):
    for neighbor, weight in zip(de_bruijn_graph[node], weighted_adjacency_matrix[node_index]):
        neighbor_index = list(de_bruijn_graph.keys()).index(neighbor)
        if labels[node_index] != labels[neighbor_index]:
            edges_cut_off += weight
            if weight > 0:
                cut_off_edges_count += 1
                if weight > 0:
                    cut_off_edges_with_weight_gt_zero += 1

cut_off_edges = []

for node_index, node in enumerate(de_bruijn_graph.keys()):
    for neighbor in de_bruijn_graph[node]:
        neighbor_index = list(de_bruijn_graph.keys()).index(neighbor)
        if labels[node_index] != labels[neighbor_index]:
            cut_off_edges.append((node, neighbor))

# Print the number of edges cut off due to partitioning
print("Number of edges cut off with weight > 0:", cut_off_edges_count)
print("Number of edges cut off with weight > 0:", cut_off_edges_with_weight_gt_zero)

# # Print the cut-off edges due to partitioning
# print("\nCut-off edges due to partitioning:")
# for edge in cut_off_edges:
#     print(edge)

# Function to partition the graph into subgraphs with 20 nodes or fewer
def partition_graph(graph, max_nodes_per_partition):
    subgraphs = []
    nodes = list(graph.nodes)
    while nodes:
        subgraph_nodes = nodes[:max_nodes_per_partition]
        subgraph = graph.subgraph(subgraph_nodes).copy()
        subgraphs.append(subgraph)
        nodes = nodes[max_nodes_per_partition:]
    return subgraphs

# Partition the graph into subgraphs with 20 nodes or fewer
subgraphs = partition_graph(G, 20)

# Print the scaffolds for each subgraph obtained
def scaffold_subgraph(subgraph):
    scaffolds = []
    for component in nx.strongly_connected_components(subgraph):
        component_subgraph = subgraph.subgraph(component)
        if len(component_subgraph.nodes) > 1:
            path = list(nx.dfs_edges(component_subgraph, source=list(component_subgraph.nodes)[0]))
            scaffold = ''.join([edge[0][0] for edge in path]) + path[-1][1]
            scaffolds.append(scaffold)
    return scaffolds

# Print scaffolds for each subgraph
for i, subgraph in enumerate(subgraphs):
    scaffolds = scaffold_subgraph(subgraph)
    print(f"\nSubgraph {i+1} scaffolds:")
    for scaffold in scaffolds:
        print(scaffold)

# Visualize each subgraph and count the number of nodes
for i, subgraph in enumerate(subgraphs):
    print(f"\nSubgraph {i+1} has {len(subgraph.nodes)} nodes")
    plt.figure()
    nx.draw(subgraph, with_labels=True, node_size=700, node_color="skyblue", font_size=10, font_weight="bold")
    plt.title(f"Subgraph {i+1}")
    plt.show()
