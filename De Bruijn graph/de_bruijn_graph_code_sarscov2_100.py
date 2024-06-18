#MAKE THE DE BRUIJN GRAPH
import networkx as nx
import matplotlib.pyplot as plt

def de_bruijn(fragments, k):
    edges = []
    nodes = set()
    kmer_count = {}  # Dictionary to store the count of each k-mer
    for fragment in fragments:
        for i in range(len(fragment) - k + 1):
            kmer = fragment[i:i+k]
            prefix = kmer[:-1]
            suffix = kmer[1:]
            edges.append((prefix, suffix))
            nodes.add(prefix)
            nodes.add(suffix)
            # Increment the count for the current k-mer
            kmer_count[kmer] = kmer_count.get(kmer, 0) + 1

    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    # Add edges with weights (frequency of k-mer)
    for edge in edges:
        # Calculate the k-mer represented by this edge
        kmer = edge[0] + edge[1][-1]
        weight = kmer_count[kmer]
        G.add_edge(edge[0], edge[1], weight=weight)

    return G

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

# Read sequences from a FASTA file
file_path = r"/content/sars_cov2_100.fasta"
sequences = read_fasta_file(file_path)

# Parameters
k = 4  # Adjust the value of k as needed

# Generate the De Bruijn graph
G = de_bruijn(sequences, k)

# Visualize the De Bruijn graph
visualize_de_bruijn(G)

# Count the number of nodes
num_nodes = G.number_of_nodes()
print("Number of nodes in the graph:", num_nodes)

# Count the total number of edges
num_edges = G.number_of_edges()
print("Number of edges in the graph:", num_edges)
