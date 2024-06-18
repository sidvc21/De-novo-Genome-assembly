import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh

# SPECTRAL PARTITIONING
def spectral_partitioning(G):
    A = nx.adjacency_matrix(G).astype(float)
    D = np.diag(np.array(A.sum(axis=1)).flatten())
    L = D - A
    eigenvalues, eigenvectors = eigsh(L, k=2, which='SM')
    second_eigenvector = eigenvectors[:, 1]
    partition1 = [node for i, node in enumerate(G.nodes()) if second_eigenvector[i] >= 0]
    partition2 = [node for i, node in enumerate(G.nodes()) if second_eigenvector[i] < 0]
    return partition1, partition2

# RECURSIVE PARTITIONING
def partition_graph_recursive(G, max_nodes=30, current_subgraph_name="Subgraph", partition_level=0):
    if G.number_of_nodes() <= max_nodes:
        return [(current_subgraph_name, G)]
    partition1, partition2 = spectral_partitioning(G)
    partition_level += 1
    subgraph1 = G.subgraph(partition1).copy()
    subgraph2 = G.subgraph(partition2).copy()
    subgraph_name_1 = f"{current_subgraph_name} Partition-{partition_level} 1"
    subgraph_name_2 = f"{current_subgraph_name} Partition-{partition_level} 2"
    subgraphs = []
    subgraphs.extend(partition_graph_recursive(subgraph1, max_nodes, subgraph_name_1, partition_level))
    subgraphs.extend(partition_graph_recursive(subgraph2, max_nodes, subgraph_name_2, partition_level))
    return subgraphs

# PRINTING SUBGRAPHS
def print_subgraphs(subgraphs):
    for idx, (subgraph_name, subgraph) in enumerate(subgraphs):
        num_nodes = subgraph.number_of_nodes()
        print(f"Subgraph {idx + 1} - Total Nodes: {num_nodes}")
        print("Node 1\tNode 2\tEdge Weight")
        for u, v, data in subgraph.edges(data=True):
            weight = data['weight'] if 'weight' in data else 1  # Default weight if not present
            print(f"{u}\t{v}\t{weight}")
        print("\n")

        # Visualize subgraph
        plt.figure(figsize=(8, 6))
        pos = nx.spring_layout(subgraph)
        edge_labels = {(u, v): d['weight'] if 'weight' in d else 1 for u, v, d in subgraph.edges(data=True)}  # Default weight if not present
        nx.draw(subgraph, pos, with_labels=True, node_size=1200, node_color='lightblue', font_size=12, font_weight='bold', edge_color='gray', arrowsize=20)
        nx.draw_networkx_edge_labels(subgraph, pos, edge_labels=edge_labels)
        plt.title(f'Subgraph {idx + 1} Visualization')
        plt.show()

# MAIN EXECUTION
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

# De Bruijn Graph Class
class DeBruijnGraph:
    def __init__(self, k):
        self.k = k
        self.graph = nx.DiGraph()

    def add_sequence(self, sequence):
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i:i+self.k]
            prefix = kmer[:-1]
            suffix = kmer[1:]
            if self.graph.has_edge(prefix, suffix):
                self.graph[prefix][suffix]['weight'] += 1
            else:
                self.graph.add_edge(prefix, suffix, weight=1)
                self.graph.nodes[prefix]['sequence'] = prefix
                self.graph.nodes[suffix]['sequence'] = suffix

    def get_graph(self):
        return self.graph

# Function to build De Bruijn graph from subgraph
def build_debruijn_graph_from_subgraph(subgraph, k):
    graph = DeBruijnGraph(k)
    for u, v, data in subgraph.edges(data=True):
        if 'sequence' in subgraph.nodes[u] and 'sequence' in subgraph.nodes[v]:
            kmer = subgraph.nodes[u]['sequence'] + subgraph.nodes[v]['sequence'][-1]
            graph.add_sequence(kmer)
    return graph.get_graph()

# Function to traverse De Bruijn graph and get contigs
def traverse_debruijn_graph(G):
    contigs = []
    for node in G.nodes:
        if G.in_degree(node) == 0 or G.out_degree(node) == 0:
            contig = node
            current = node
            while G.out_degree(current) > 0:
                next_node = list(G.successors(current))[0]
                contig += next_node[-1]
                current = next_node
            contigs.append(contig)
    return contigs

# Scaffolder Class
class Scaffolder:
    def __init__(self, contigs):
        self.contigs = contigs
        self.scaffolds = []

    def scaffold_contigs(self):
        self.scaffolds.append(self.contigs[0])
        for contig in self.contigs[1:]:
            merged = False
            for scaffold in self.scaffolds:
                overlap = self.find_overlap(scaffold, contig)
                if overlap:
                    self.extend_scaffold(scaffold, contig, overlap)
                    merged = True
                    break
            if not merged:
                self.scaffolds.append(contig)

    def find_overlap(self, scaffold, contig):
        k = min(len(scaffold), len(contig))
        for i in range(k, 0, -1):
            if scaffold.endswith(contig[:i]):
                return i
        return 0

    def extend_scaffold(self, scaffold, contig, overlap):
        self.scaffolds[self.scaffolds.index(scaffold)] += contig[overlap:]

    def print_and_store_scaffolds(self, output_file):
        with open(output_file, "w") as f:
            contig_num = 1
            for contig in self.contigs:
                f.write(f"Contig {contig_num}:\n")
                f.write(f"{contig}\n\n")
                contig_num += 1
            scaffold_num = 1
            for scaffold in self.scaffolds:
                f.write(f"Scaffold {scaffold_num}:\n")
                f.write(f"{scaffold}\n\n")

def de_bruijn(sequences, k):
    dbg = DeBruijnGraph(k)
    for seq in sequences:
        dbg.add_sequence(seq)
    return dbg.get_graph()

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

# Processing Subgraphs
for idx, (subgraph_name, subgraph) in enumerate(subgraphs):
    subgraph_k = k - 1  # Adjust k-mer size for subgraph De Bruijn graph
    subgraph_debruijn = build_debruijn_graph_from_subgraph(subgraph, subgraph_k)
    contigs = traverse_debruijn_graph(subgraph_debruijn)
    scaffolder = Scaffolder(contigs)
    scaffolder.scaffold_contigs()
    output_file = f"subgraph_{idx+1}_scaffolds.txt"
    scaffolder.print_and_store_scaffolds(output_file)
    print(f"Scaffolds for {subgraph_name} saved to {output_file}")
