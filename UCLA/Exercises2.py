def string_composition(k, text):
    """
    Generate all k-mers from the input text.

    Parameters:
    k (int): Length of each k-mer.
    text (str): Input string.

    Returns:
    list: List of k-mers.
    """
    # Create a list of all k-mers in the text
    kmers = [text[i:i+k] for i in range(len(text) - k + 1)]
    return kmers

# Define the file path to the dataset
file_path = 'dataset_30153_3.txt'

# Read the input dataset
with open(file_path, 'r') as file:
    # Extract the input from the file
    k = int(file.readline().strip())
    text = file.readline().strip()

# Compute the k-mer composition
result = string_composition(k, text)

# Format the output as a space-separated string
output = " ".join(result)

# Print the output
print(output)

# Optionally, save the result to a file
with open('output.txt', 'w') as output_file:
    output_file.write(output)


############

# Function to reconstruct the string from the genome path
def reconstruct_string(genome_path):
    """
    Reconstructs a string from a genome path of overlapping k-mers.

    Parameters:
        genome_path (list of str): List of overlapping k-mers.

    Returns:
        str: Reconstructed string.
    """
    # Initialize the resulting string with the first k-mer
    result = genome_path[0]
    
    # Iterate through the rest of the k-mers
    for kmer in genome_path[1:]:
        # Append the last character of the current k-mer
        result += kmer[-1]
    
    return result

# Load the dataset from the provided file
file_path = 'dataset_30182_3.txt'

# Read the file and split into genome path k-mers
with open(file_path, 'r') as file:
    genome_path = file.read().strip().split()

# Reconstruct the string from the genome path
output_string = reconstruct_string(genome_path)

# Save the output to a file
output_file = 'reconstructed_genome_string.txt'
with open(output_file, 'w') as file:
    file.write(output_string)

# Print confirmation and the path to the output file
print(f"The reconstructed genome string has been saved to: {output_file}")
##################


# Function to build the overlap graph for a collection of k-mers
def build_overlap_graph(patterns):
    """
    Constructs the overlap graph for a collection of k-mers.

    Parameters:
        patterns (list of str): A list of k-mers.

    Returns:
        dict: A dictionary representing the adjacency list of the overlap graph.
    """
    adjacency_list = {}

    # Create mappings for prefixes and suffixes
    prefix_dict = {}
    suffix_dict = {}

    for pattern in patterns:
        prefix = pattern[:-1]
        suffix = pattern[1:]
        prefix_dict[pattern] = prefix
        suffix_dict[pattern] = suffix

    # Construct the overlap graph
    for a in patterns:
        adjacency = []
        for b in patterns:
            if a != b and suffix_dict[a] == prefix_dict[b]:
                adjacency.append(b)
        if adjacency:
            adjacency_list[a] = adjacency

    return adjacency_list

# Load the dataset from the provided file
file_path = 'dataset_30182_10.txt'

# Read the file and split into genome path k-mers
with open(file_path, 'r') as file:
    genome_path = file.read().strip().split()

# Build the overlap graph
overlap_graph = build_overlap_graph(genome_path)

# Save the adjacency list to a file
output_file = 'overlap_graph_output.txt'
with open(output_file, 'w') as file:
    for node in sorted(overlap_graph.keys()):
        edges = ' '.join(overlap_graph[node])
        file.write(f"{node}: {edges}\n")

print(f"The overlap graph has been saved to: {output_file}")

#########################

from itertools import product
from collections import defaultdict

def generate_binary_strings(k):
    """
    Generate all binary strings of length k.
    """
    return [''.join(seq) for seq in product('01', repeat=k)]

def build_de_bruijn_graph(k):
    """
    Build the de Bruijn graph for binary strings of length k-1.
    """
    graph = defaultdict(list)
    for binary in generate_binary_strings(k):
        prefix = binary[:-1]
        suffix = binary[1:]
        graph[prefix].append(suffix)
    return graph

def eulerian_cycle(graph):
    """
    Find an Eulerian cycle in the de Bruijn graph.
    """
    stack = []
    cycle = []
    current = next(iter(graph))  # Start at any node
    while stack or graph[current]:
        if not graph[current]:
            cycle.append(current)
            current = stack.pop()
        else:
            stack.append(current)
            next_node = graph[current].pop()
            current = next_node
    cycle.append(current)
    return cycle[::-1]  # Reverse to get the correct order

def construct_universal_string(k):
    """
    Construct a k-universal string using the de Bruijn sequence.
    """
    graph = build_de_bruijn_graph(k)
    cycle = eulerian_cycle(graph)
    # Generate the universal string by concatenating nodes from the Eulerian cycle
    universal_string = ''.join(node[0] for node in cycle[:-1]) + cycle[-1]
    return universal_string

# Generate the 4-universal string
k = 4
universal_string = construct_universal_string(k)

# Output the result
print(f"The 4-universal string is: {universal_string}")

# Verify that the string contains all k-mers
def verify_universal_string(universal_string, k):
    expected = set(generate_binary_strings(k))
    observed = {universal_string[i:i+k] for i in range(len(universal_string) - k + 1)}
    return expected == observed

# Verification step
is_correct = verify_universal_string(universal_string, k)
print(f"Does the universal string contain all 4-mers? {is_correct}")
###################



from collections import defaultdict

def de_bruijn_graph(k, text):
    """
    Constructs the de Bruijn graph from a given string and integer k.
    
    Parameters:
        k (int): Length of k-mers.
        text (str): Input string.
    
    Returns:
        dict: Adjacency list of the de Bruijn graph.
    """
    adjacency_list = defaultdict(list)
    
    # Generate all k-mers from the text
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        prefix = kmer[:-1]
        suffix = kmer[1:]
        adjacency_list[prefix].append(suffix)
    
    # Format the adjacency list into a dictionary
    formatted_adjacency = {node: sorted(targets) for node, targets in adjacency_list.items()}
    
    return formatted_adjacency

# Example input
k = 4
text = "AAGATTCTCTAAGA"

# Generate the de Bruijn graph
graph = de_bruijn_graph(k, text)

# Print the adjacency list in the required format
for node, edges in graph.items():
    print(f"{node}: {' '.join(edges)}")
###########################


from collections import defaultdict

def de_bruijn_graph(k, text):
    """
    Constructs the De Bruijn graph from a given string and integer k.
    
    Parameters:
        k (int): Length of k-mers.
        text (str): Input string.
    
    Returns:
        dict: Adjacency list of the De Bruijn graph.
    """
    adjacency_list = defaultdict(list)
    
    # Generate all k-mers from the text
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        prefix = kmer[:-1]
        suffix = kmer[1:]
        adjacency_list[prefix].append(suffix)
    
    # Format the adjacency list into a dictionary with sorted values for consistent output
    formatted_adjacency = {node: sorted(targets) for node, targets in adjacency_list.items()}
    
    return formatted_adjacency

# Input and output file paths
input_file = "dataset_30183_6.txt"
output_file = "de_bruijn_graph_output.txt"

# Read input
with open(input_file, 'r') as file:
    lines = file.readlines()
    k = int(lines[0].strip())
    text = lines[1].strip()

# Generate the De Bruijn graph
graph = de_bruijn_graph(k, text)

# Write output
with open(output_file, 'w') as file:
    for node, edges in graph.items():
        file.write(f"{node}: {' '.join(edges)}\n")

print(f"The De Bruijn graph has been saved to {output_file}")
#########################

from collections import defaultdict

def de_bruijn_from_kmers(patterns):
    """
    Constructs the De Bruijn graph from a collection of k-mers.
    
    Parameters:
        patterns (list of str): Collection of k-mers.
    
    Returns:
        dict: Adjacency list of the De Bruijn graph.
    """
    adjacency_list = defaultdict(list)
    
    # Construct the adjacency list
    for kmer in patterns:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        adjacency_list[prefix].append(suffix)
    
    # Format the adjacency list for consistent output
    formatted_adjacency = {node: sorted(targets) for node, targets in adjacency_list.items()}
    return formatted_adjacency

# File paths
input_file = "dataset_30184_8.txt"  # Replace with the actual input file name
output_file = "de_bruijn_graph_from_kmers_output.txt"  # Replace with the desired output file name

# Load k-mers from the input file
with open(input_file, 'r') as file:
    patterns = file.read().strip().split()

# Generate the De Bruijn graph
graph = de_bruijn_from_kmers(patterns)

# Save the adjacency list to the output file
with open(output_file, 'w') as file:
    for node, edges in graph.items():
        file.write(f"{node}: {' '.join(edges)}\n")

print(f"The De Bruijn graph has been saved to: {output_file}")
############## Test
# Adjacency list
graph = {
    1: [2, 3, 5],
    2: [4, 5],
    3: [1, 2, 5],
    4: [1, 3],
    5: [2, 4]
}

def is_hamiltonian_cycle(graph, cycle):
    """
    Check if the given cycle is a Hamiltonian cycle.
    
    Parameters:
        graph (dict): The adjacency list of the graph.
        cycle (list): The cycle to check.
    
    Returns:
        bool: True if the cycle is Hamiltonian, False otherwise.
    """
    # Check if the cycle visits every node exactly once
    nodes = set(graph.keys())
    if len(cycle) != len(nodes) + 1 or cycle[0] != cycle[-1]:
        return False

    # Check if every consecutive pair of nodes in the cycle has an edge
    for i in range(len(cycle) - 1):
        if cycle[i + 1] not in graph[cycle[i]]:
            return False

    return True

# Given cycles
cycles = [
    [1, 3, 2, 5, 4, 1],
    [1, 5, 4, 3, 2, 1],
    [1, 2, 4, 3, 5, 1],
    [2, 4, 3, 5, 2]
]

# Check each cycle
for cycle in cycles:
    print(f"Cycle {cycle} is Hamiltonian: {is_hamiltonian_cycle(graph, cycle)}")


######
from itertools import product

def is_3_universal(string, k=3):
    """
    Check if the given string is k-universal.
    
    Parameters:
        string (str): The binary string to check.
        k (int): The length of substrings to check.
    
    Returns:
        bool: True if the string is k-universal, False otherwise.
    """
    substrings = {string[i:i + k] for i in range(len(string) - k + 1)}
    all_kmers = {''.join(p) for p in product('01', repeat=k)}
    return substrings == all_kmers

# Given strings
strings = [
    "1110001011",
    "0011100100",
    "0101010100",
    "1100011011",
    "1101000111",
    "1000101110"
]

# Check each string
for s in strings:
    print(f"String {s} is 3-universal: {is_3_universal(s)}")
#####


from collections import defaultdict

def calculate_indegree(kmers, target):
    """
    Calculate the indegree of a given node in the De Bruijn graph.
    
    Parameters:
        kmers (list): List of k-mers.
        target (str): The target node to calculate the indegree for.
    
    Returns:
        int: The indegree of the target node.
    """
    indegree = 0
    for kmer in kmers:
        suffix = kmer[1:]  # Get the suffix of the k-mer
        if suffix == target:
            indegree += 1
    return indegree

kmers = [
    "GCGA", "CAAG", "AAGA", "GCCG", "ACAA", "AGTA", "TAGG", "AGTA",
    "ACGT", "AGCC", "TTCG", "AGTT", "AGTA", "CGTA", "GCGC", "GCGA",
    "GGTC", "GCAT", "AAGC", "TAGA", "ACAG", "TAGA", "TCCT", "CCCC",
    "GCGC", "ATCC", "AGTA", "AAGA", "GCGA", "CGTA"
]

# Target node
target = "GTA"

# Calculate indegree
indegree = calculate_indegree(kmers, target)
print(f"The indegree of {target} is: {indegree}")

#################### Part 2 Exercise

#1

from collections import defaultdict, deque

def find_eulerian_path(graph):
    """Finds an Eulerian path in a directed graph."""
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    for node in graph:
        out_degree[node] += len(graph[node])
        for neighbor in graph[node]:
            in_degree[neighbor] += 1
    
    # Find start and end nodes
    start, end = None, None
    for node in set(in_degree.keys()).union(out_degree.keys()):
        if out_degree[node] - in_degree[node] == 1:
            start = node
        elif in_degree[node] - out_degree[node] == 1:
            end = node
    
    if not start:
        start = next(iter(graph))  # Arbitrary start for Eulerian cycle

    # Hierholzer's Algorithm
    stack = [start]
    path = []
    while stack:
        node = stack[-1]
        if graph[node]:
            next_node = graph[node].pop()
            stack.append(next_node)
        else:
            path.append(stack.pop())
    
    return path[::-1]

def assemble_sequence(kmers):
    """Assembles a sequence from a list of k-mers."""
    # Build the de Bruijn graph
    graph = defaultdict(list)
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph[prefix].append(suffix)
    
    # Find Eulerian path
    path = find_eulerian_path(graph)
    
    # Reconstruct sequence from path
    sequence = path[0]
    for node in path[1:]:
        sequence += node[-1]
    
    return sequence

# Input 4-mers
kmers = [
    "AAAT", "AATG", "ACCC", "ACGC", "ATAC", "ATCA", "ATGC",
    "CAAA", "CACC", "CATA", "CATC", "CCAG", "CCCA", "CGCT",
    "CTCA", "GCAT", "GCTC", "TACG", "TCAC", "TCAT", "TGCA"
]


# Assemble the sequence
result = assemble_sequence(kmers)
print("Assembled Sequence:", result)




#2 
from collections import defaultdict

def calculate_minimum_edges(adj_list):
    # Step 1: Calculate in-degree and out-degree
    out_degree = defaultdict(int)
    in_degree = defaultdict(int)

    for node, neighbors in adj_list.items():
        out_degree[node] += len(neighbors)
        for neighbor in neighbors:
            in_degree[neighbor] += 1

    # Step 2: Calculate imbalance for each node
    imbalance = defaultdict(int)
    for node in set(out_degree.keys()).union(in_degree.keys()):
        imbalance[node] = out_degree[node] - in_degree[node]

    # Step 3: Count total positive imbalance
    total_positive_imbalance = sum(max(0, imbalance[node]) for node in imbalance)

    return total_positive_imbalance

# Adjacency list of the graph
adj_list = {
    1: [2, 3, 5],
    2: [1, 4],
    3: [2, 5],
    4: [1, 2, 5],
    5: [3]
}

# Calculate the minimum number of edges to balance the graph
min_edges = calculate_minimum_edges(adj_list)
print("Minimum number of edges to add:", min_edges)



# 3
# Correcting the logic for reconstructing the string from the (3,1)-mer composition

from collections import defaultdict, deque

def build_graph(kmers):
    """Builds the de Bruijn graph from the (3,1)-mers."""
    graph = defaultdict(list)
    for prefix, suffix in kmers:
        graph[prefix].append(suffix)
    return graph

def find_eulerian_path(graph):
    """Finds an Eulerian path in the graph using Hierholzer's algorithm."""
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)

    for node in graph:
        out_degree[node] += len(graph[node])
        for neighbor in graph[node]:
            in_degree[neighbor] += 1

    # Find start and end nodes
    start, end = None, None
    for node in set(in_degree.keys()).union(out_degree.keys()):
        if out_degree[node] - in_degree[node] == 1:
            start = node
        elif in_degree[node] - out_degree[node] == 1:
            end = node

    if not start:
        start = next(iter(graph))  # Arbitrary start if Eulerian cycle

    # Hierholzer's algorithm
    stack = [start]
    path = []
    while stack:
        node = stack[-1]
        if graph[node]:
            next_node = graph[node].pop()
            stack.append(next_node)
        else:
            path.append(stack.pop())

    return path[::-1]

def reconstruct_string(kmers):
    """Reconstructs the string from the (3,1)-mer composition."""
    # Build the de Bruijn graph
    graph = build_graph(kmers)

    # Find the Eulerian path
    path = find_eulerian_path(graph)

    # Reconstruct the string
    sequence = path[0]
    for node in path[1:]:
        sequence += node[-1]

    return sequence

# Input (3,1)-mers
kmers = [
    ("ACC", "ATA"), ("ACT", "ATT"), ("ATA", "TGA"), ("ATT", "TGA"),
    ("CAC", "GAT"), ("CCG", "TAC"), ("CGA", "ACT"), ("CTG", "AGC"),
    ("CTG", "TTC"), ("GAA", "CTT"), ("GAT", "CTG"), ("GAT", "CTG"),
    ("TAC", "GAT"), ("TCT", "AAG"), ("TGA", "GCT"), ("TGA", "TCT"),
    ("TTC", "GAA")
]

# Reconstruct the string
result = reconstruct_string(kmers)
print("Reconstructed String:", result)

############################
def find_eulerian_cycle(graph):
    """
    Finds an Eulerian cycle in a directed graph using Hierholzer's algorithm.
    
    :param graph: Dictionary where keys are nodes and values are lists of neighbors
    :return: List of nodes in Eulerian cycle
    """
    # Step 1: Create a stack to hold the current path and a list for the final cycle
    stack = []
    cycle = []

    # Step 2: Start with any node (e.g., the first key in the graph)
    current_node = next(iter(graph))
    stack.append(current_node)

    while stack:
        if graph[current_node]:
            # If the current node has neighbors, push it to the stack and follow an edge
            stack.append(current_node)
            next_node = graph[current_node].pop()
            current_node = next_node
        else:
            # If the current node has no neighbors, add it to the cycle and backtrack
            cycle.append(current_node)
            current_node = stack.pop()

    # Reverse the cycle to get the correct order
    return cycle[::-1]

# Reading the adjacency list from the input
def parse_input(input_file):
    """
    Parses the input file into a graph represented as an adjacency list.
    :param input_file: Path to the input file
    :return: Dictionary representing the adjacency list
    """
    graph = {}
    with open(input_file, 'r') as file:
        for line in file:
            parts = line.strip().split(': ')
            node = int(parts[0])
            neighbors = list(map(int, parts[1].split()))
            graph[node] = neighbors
    return graph

# Example usage
if __name__ == "__main__":
    # Replace 'dataset_30187_2.txt' with your actual file name
    input_file = "dataset_30187_2.txt"
    graph = parse_input(input_file)
    eulerian_cycle = find_eulerian_cycle(graph)
    print(" ".join(map(str, eulerian_cycle)))
#######################

from collections import defaultdict, deque

def find_eulerian_path(graph):
    # Helper function to find the starting point
    def find_start_node(graph):
        in_degree = defaultdict(int)
        out_degree = defaultdict(int)
        for node, neighbors in graph.items():
            out_degree[node] += len(neighbors)
            for neighbor in neighbors:
                in_degree[neighbor] += 1
        
        start, end = None, None
        for node in set(in_degree.keys()).union(out_degree.keys()):
            if out_degree[node] - in_degree[node] == 1:
                start = node
            elif in_degree[node] - out_degree[node] == 1:
                end = node
        return start if start else list(graph.keys())[0]
    
    # Hierholzer's Algorithm to find the Eulerian Path
    def hierholzer(graph, start):
        stack = [start]
        path = []
        current_path = deque()
        while stack:
            while graph[stack[-1]]:
                neighbor = graph[stack[-1]].pop()
                stack.append(neighbor)
            current_path.append(stack.pop())
        return list(current_path)[::-1]

    # Copy graph to avoid modifying the input
    graph_copy = {node: deque(neighbors) for node, neighbors in graph.items()}
    start_node = find_start_node(graph_copy)
    return hierholzer(graph_copy, start_node)

# Parse the input graph
def parse_input(input_str):
    graph = defaultdict(list)
    for line in input_str.strip().split("\n"):
        node, neighbors = line.split(": ")
        graph[int(node)] = list(map(int, neighbors.split()))
    return graph

# Example Input
input_str = """0: 3
1: 6
3: 1
4: 1
5: 4
6: 5 7
7: 9
8: 6
9: 8"""

graph = parse_input(input_str)
eulerian_path = find_eulerian_path(graph)
print(" ".join(map(str, eulerian_path)))
#####################

from collections import defaultdict

def build_de_bruijn_graph(kmers):
    """Builds the de Bruijn graph from the k-mers"""
    graph = defaultdict(list)
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    
    # Iterate through each k-mer to build the graph
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph[prefix].append(suffix)
        out_degree[prefix] += 1
        in_degree[suffix] += 1

    return graph, in_degree, out_degree

def find_eulerian_path(graph, in_degree, out_degree):
    """Finds the Eulerian path using Hierholzer's algorithm"""
    start_node = None
    for node in graph:
        if out_degree[node] > in_degree[node]:
            start_node = node
            break
    
    if not start_node:
        start_node = next(iter(graph))
    
    path = []
    stack = [start_node]
    
    while stack:
        node = stack[-1]
        if graph[node]:
            next_node = graph[node].pop()
            stack.append(next_node)
        else:
            path.append(stack.pop())
    
    return path[::-1]

def path_to_genome(path, k):
    """Reconstructs the genome from the Eulerian path"""
    genome = path[0]
    for node in path[1:]:
        genome += node[-1]
    return genome

def string_reconstruction(kmers):
    k = len(kmers[0])
    graph, in_degree, out_degree = build_de_bruijn_graph(kmers)
    eulerian_path = find_eulerian_path(graph, in_degree, out_degree)
    genome = path_to_genome(eulerian_path, k)
    return genome

# Read the input from the user or a file
with open('dataset_30187_7.txt', 'r') as file:
    data = file.readlines()

k = int(data[0].strip())
kmers = data[1].strip().split()

# Solve the problem
result = string_reconstruction(kmers)

# Output the result
print(result)

######################################
from collections import defaultdict
import itertools

def generate_unique_binary_kmers(k):
    """Generate all unique binary k-mers."""
    return [''.join(seq) for seq in itertools.product('01', repeat=k)]

def debruijn_graph(kmers):
    """Construct the de Bruijn graph from k-mers."""
    graph = defaultdict(list)
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph[prefix].append(suffix)
    return graph

def eulerian_cycle(graph):
    """Find an Eulerian cycle in a directed graph."""
    stack = []
    path = []
    # Start with any node that has outgoing edges
    start_node = next(iter(graph))
    stack.append(start_node)

    while stack:
        current = stack[-1]
        if graph[current]:
            stack.append(graph[current].pop())
        else:
            path.append(stack.pop())
    return path[::-1]  # Reverse to get the correct order

def glue_nodes_correctly(cycle, k):
    """Glue nodes to form a circular string by skipping the last (k-1) nodes."""
    circular_string = cycle[0]  # Start with the first node
    for node in cycle[1:]:
        circular_string += node[-1]
    # Skip the last (k-1) characters to make it circular
    return circular_string[:-(k-1)]

def universal_string(k):
    """Generate a k-universal circular string."""
    str_kmers = generate_unique_binary_kmers(k)
    dbjn_graph = debruijn_graph(str_kmers)
    cycle = eulerian_cycle(dbjn_graph)
    glued = glue_nodes_correctly(cycle, k)

    # Print or return the result
    return glued

# Example usage
k = 8
result = universal_string(k)
print(result)

#######################

def generate_k_d_mer_composition(sequence, k, d):
    """Generate the (k,d)-mer composition of a sequence."""
    composition = []
    for i in range(len(sequence) - k - d + 1):
        prefix = sequence[i:i+k]
        suffix = sequence[i+k+d:i+2*k+d]
        composition.append(f"({prefix}|{suffix})")
    return sorted(composition)  # Return in lexicographic order

# Input sequence and parameters
sequence = "TAATGCCATGGGATGTT"
k = 3
d = 2

# Generate (k,d)-mer composition
composition = generate_k_d_mer_composition(sequence, k, d)

# Join the result into a single line
result = " ".join(composition)
result
#####################






