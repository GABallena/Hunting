def HammingDistance(p, q):
    """
    Computes the Hamming distance between two strings of equal length.
    """
    return sum(1 for pi, qi in zip(p, q) if pi != qi)

def ApproximatePatternMatching(Pattern, Text, d):
    """
    Finds all starting positions where Pattern appears as a substring of Text with at most d mismatches.

    Parameters:
        Pattern (str): The pattern to search for.
        Text (str): The string to search within.
        d (int): The maximum number of mismatches allowed.

    Returns:
        list: A list of starting positions (0-based) where Pattern appears with at most d mismatches.
    """
    positions = []  # List to store the starting positions
    k = len(Pattern)

    for i in range(len(Text) - k + 1):
        substring = Text[i:i + k]
        if HammingDistance(Pattern, substring) <= d:
            positions.append(i)

    return positions

# Example of how to load and process input
def process_input_file(file_path):
    """
    Reads the dataset and extracts the Pattern, Text, and d values.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
        Pattern = lines[0].strip()
        Text = lines[1].strip()
        d = int(lines[2].strip())
    return Pattern, Text, d

# File path (replace with the actual dataset path)
file_path = 'dataset_30278_6.txt'

# Load the data
Pattern, Text, d = process_input_file(file_path)

# Solve the problem
positions = ApproximatePatternMatching(Pattern, Text, d)

# Print the positions
print(" ".join(map(str, positions)))

#####
def HammingDistance(p, q):
    """
    Computes the Hamming distance between two strings of equal length.
    """
    return sum(1 for pi, qi in zip(p, q) if pi != qi)

def ApproximatePatternMatching(Pattern, Text, d):
    """
    Finds all starting positions where Pattern appears as a substring of Text with at most d mismatches.

    Parameters:
        Pattern (str): The pattern to search for.
        Text (str): The string to search within.
        d (int): The maximum number of mismatches allowed.

    Returns:
        list: A list of starting positions (0-based) where Pattern appears with at most d mismatches.
    """
    positions = []  # List to store the starting positions
    k = len(Pattern)

    for i in range(len(Text) - k + 1):
        substring = Text[i:i + k]
        if HammingDistance(Pattern, substring) <= d:
            positions.append(i)

    return positions

# Input values
Text = "AACAAGCTGATAAACATTTAAAGAG"
Pattern = "AAAA"
d = 2

# Find approximate matches
positions = ApproximatePatternMatching(Pattern, Text, d)

# Count the number of items in the positions array
count = len(positions)
print("Count with at most 2 mismatches:", count)

###################

# Re-confirm with Python
Text = "AACAAGCTGATAAACATTTAAAGAG"
Pattern = "AAAAA"
d = 2

# Function Definitions
def HammingDistance(p, q):
    return sum(1 for pi, qi in zip(p, q) if pi != qi)

def Countd(Text, Pattern, d):
    count = 0
    k = len(Pattern)
    for i in range(len(Text) - k + 1):
        substring = Text[i:i + k]
        if HammingDistance(Pattern, substring) <= d:
            count += 1
    return count

# Compute Count2
count = Countd(Text, Pattern, d)
print("Count2:", count)
#########################################################


from itertools import product

def HammingDistance(p, q):
    """
    Computes the Hamming distance between two strings of equal length.
    """
    return sum(1 for pi, qi in zip(p, q) if pi != qi)

def Neighbors(Pattern, d):
    """
    Generates all k-mers within d mismatches of Pattern.
    """
    nucleotides = "ACGT"
    if d == 0:
        return {Pattern}
    if len(Pattern) == 0:
        return {""}

    first_symbol = Pattern[0]
    suffix = Pattern[1:]
    suffix_neighbors = Neighbors(suffix, d)
    neighbors = set()

    for text in suffix_neighbors:
        if HammingDistance(suffix, text) < d:
            for nucleotide in nucleotides:
                neighbors.add(nucleotide + text)
        else:
            neighbors.add(first_symbol + text)

    return neighbors

def FrequentWordsWithMismatches(Text, k, d):
    """
    Finds the most frequent k-mers with up to d mismatches in Text.
    """
    patterns = {}
    n = len(Text)

    for i in range(n - k + 1):
        current_pattern = Text[i:i + k]
        neighbors = Neighbors(current_pattern, d)
        for neighbor in neighbors:
            patterns[neighbor] = patterns.get(neighbor, 0) + 1

    max_count = max(patterns.values())
    most_frequent = [pattern for pattern, count in patterns.items() if count == max_count]

    return most_frequent

def load_dataset(file_path):
    """
    Loads the dataset and extracts the input values.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
        Text = lines[0].strip()
        k, d = map(int, lines[1].strip().split())
    return Text, k, d

# Load the dataset
file_path = "dataset_30278_9.txt"  # Replace with your dataset's actual path
Text, k, d = load_dataset(file_path)

# Compute the result
result = FrequentWordsWithMismatches(Text, k, d)
print("Most frequent k-mers with mismatches:", " ".join(result))
############################

from itertools import product

def HammingDistance(p, q):
    """
    Computes the Hamming distance between two strings of equal length.
    """
    return sum(1 for pi, qi in zip(p, q) if pi != qi)

def ReverseComplement(Pattern):
    """
    Computes the reverse complement of a DNA string.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement[base] for base in reversed(Pattern))

def Neighbors(Pattern, d):
    """
    Generates all k-mers within d mismatches of Pattern.
    """
    nucleotides = "ACGT"
    if d == 0:
        return {Pattern}
    if len(Pattern) == 0:
        return {""}

    first_symbol = Pattern[0]
    suffix = Pattern[1:]
    suffix_neighbors = Neighbors(suffix, d)
    neighbors = set()

    for text in suffix_neighbors:
        if HammingDistance(suffix, text) < d:
            for nucleotide in nucleotides:
                neighbors.add(nucleotide + text)
        else:
            neighbors.add(first_symbol + text)

    return neighbors

def FrequentWordsWithMismatchesAndReverseComplements(Text, k, d):
    """
    Finds the most frequent k-mers with mismatches and reverse complements in Text.
    """
    patterns = {}
    n = len(Text)

    for i in range(n - k + 1):
        current_pattern = Text[i:i + k]
        reverse_complement = ReverseComplement(current_pattern)

        # Generate neighbors for the pattern and its reverse complement
        neighbors = Neighbors(current_pattern, d)
        neighbors_rc = Neighbors(reverse_complement, d)

        for neighbor in neighbors:
            patterns[neighbor] = patterns.get(neighbor, 0) + 1
        for neighbor_rc in neighbors_rc:
            patterns[neighbor_rc] = patterns.get(neighbor_rc, 0) + 1

    max_count = max(patterns.values())
    most_frequent = [pattern for pattern, count in patterns.items() if count == max_count]

    return most_frequent

def load_dataset(file_path):
    """
    Loads the dataset and extracts the input values.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
        Text = lines[0].strip()
        k, d = map(int, lines[1].strip().split())
    return Text, k, d

# File path (update with the actual path to your dataset)
file_path = "dataset_30278_10.txt"  # Update this with the correct file path

# Load the dataset
Text, k, d = load_dataset(file_path)

# Compute the result
result = FrequentWordsWithMismatchesAndReverseComplements(Text, k, d)
print("Most frequent k-mers with mismatches and reverse complements:", " ".join(result))
#################################

from itertools import product

def HammingDistance(p, q):
    """
    Computes the Hamming distance between two strings of equal length.
    """
    return sum(1 for pi, qi in zip(p, q) if pi != qi)

def ReverseComplement(Pattern):
    """
    Computes the reverse complement of a DNA string.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement[base] for base in reversed(Pattern))

def Neighbors(Pattern, d):
    """
    Generates all k-mers within d mismatches of Pattern.
    """
    nucleotides = "ACGT"
    if d == 0:
        return {Pattern}
    if len(Pattern) == 0:
        return {""}

    first_symbol = Pattern[0]
    suffix = Pattern[1:]
    suffix_neighbors = Neighbors(suffix, d)
    neighbors = set()

    for text in suffix_neighbors:
        if HammingDistance(suffix, text) < d:
            for nucleotide in nucleotides:
                neighbors.add(nucleotide + text)
        else:
            neighbors.add(first_symbol + text)

    return neighbors

def FrequentWordsWithMismatchesAndReverseComplements(Text, k, d):
    """
    Finds the most frequent k-mers with mismatches and reverse complements in Text.
    """
    patterns = {}
    n = len(Text)

    for i in range(n - k + 1):
        current_pattern = Text[i:i + k]
        reverse_complement = ReverseComplement(current_pattern)

        # Generate neighbors for the pattern and its reverse complement
        neighbors = Neighbors(current_pattern, d)
        neighbors_rc = Neighbors(reverse_complement, d)

        for neighbor in neighbors:
            patterns[neighbor] = patterns.get(neighbor, 0) + 1
        for neighbor_rc in neighbors_rc:
            patterns[neighbor_rc] = patterns.get(neighbor_rc, 0) + 1

    max_count = max(patterns.values())
    most_frequent = [pattern for pattern, count in patterns.items() if count == max_count]

    return most_frequent

def ComputeSkew(Genome):
    """
    Computes the skew diagram for a given genome string.

    Parameters:
        Genome (str): The DNA string.

    Returns:
        list: A list of integers representing Skew_i values from i = 0 to |Genome|.
    """
    skew = [0]
    for nucleotide in Genome:
        if nucleotide == 'G':
            skew.append(skew[-1] + 1)
        elif nucleotide == 'C':
            skew.append(skew[-1] - 1)
        else:
            skew.append(skew[-1])
    return skew

def FindOriRegion(Genome):
    """
    Finds the ori region by locating the minimum skew.
    """
    skew_values = ComputeSkew(Genome)
    min_value = min(skew_values)
    min_indices = [i for i, value in enumerate(skew_values) if value == min_value]
    return min_indices

# File path to the Salmonella_enterica genome
file_path = "Salmonella_enterica.txt"  # Replace with the correct file path

# Load the genome data (cleaning up newlines and whitespace)
with open(file_path, 'r') as file:
    Genome = file.read().replace('\n', '').strip()


# Find the ori region
ori_indices = FindOriRegion(Genome)
print("Ori region indices:", ori_indices)

# Extract a region around the minimum skew
ori_region_start = max(0, ori_indices[0] - 1000)
ori_region_end = min(len(Genome), ori_indices[0] + 1000)
ori_region = Genome[ori_region_start:ori_region_end]

# Find the DnaA box
k = 9  # Length of DnaA box
d = 1  # Allow 1 mismatch
result = FrequentWordsWithMismatchesAndReverseComplements(ori_region, k, d)
print("Candidate DnaA boxes:", " ".join(result))

#############

def HammingDistance(p, q):
    """
    Computes the Hamming distance between two strings of equal length.
    """
    return sum(1 for pi, qi in zip(p, q) if pi != qi)

def Neighbors(Pattern, d):
    """
    Generates all k-mers within d mismatches of Pattern.

    Parameters:
        Pattern (str): The DNA string.
        d (int): The maximum number of mismatches.

    Returns:
        set: A set of all strings within d mismatches of Pattern.
    """
    nucleotides = "ACGT"
    if d == 0:
        return {Pattern}
    if len(Pattern) == 0:
        return {""}

    first_symbol = Pattern[0]
    suffix = Pattern[1:]
    suffix_neighbors = Neighbors(suffix, d)
    neighbors = set()

    for text in suffix_neighbors:
        if HammingDistance(suffix, text) < d:
            for nucleotide in nucleotides:
                neighbors.add(nucleotide + text)
        else:
            neighbors.add(first_symbol + text)

    return neighbors

# Load dataset
file_path = "dataset_30282_4.txt"  # Update this path with the correct dataset location

# Read the input data
with open(file_path, 'r') as file:
    lines = file.readlines()
    Pattern = lines[0].strip()
    d = int(lines[1].strip())

# Compute the neighbors
result = Neighbors(Pattern, d)

# Print the results
print("\n".join(result))
#################################################
def HammingDist(p, q):
    """
    Calculate the Hamming distance between two sequences p and q.
    
    The Hamming distance between two strings is the number of positions at which the corresponding symbols are different.
    
    Parameters:
    p (sequence): The first sequence.
    q (sequence): The second sequence.
    
    Returns:
    int: The Hamming distance between p and q.
    
    Note:
    The lengths of the two sequences do not need to be the same. The function compares elements at corresponding positions until the end of the shorter sequence is reached.
    """
    # Use the zip function to pair up elements from both sequences and generate a list of booleans indicating whether the elements are unequal
    # Sum the booleans (True values count as 1) to calculate the total number of differences, which represents the Hamming distance
    return sum(p != q for p, q in zip(p, q))

def ImmediateNeighbor(Pattern):
    """
    Finds all immediate neighbor patterns for a given pattern.

    Parameters:
    Pattern (str): The input pattern string.

    Returns:
    list: A list containing all immediate neighbor patterns.
    """
    # Initialize the neighborhood list with the input pattern itself
    Neighborhood = [Pattern]
    # Iterate over each position in the pattern
    for i in range(len(Pattern)):
        # Iterate over all nucleotides
        for j in ["A", "C", "G", "T"]:
            # If the current nucleotide is not equal to j, replace it and add to the neighborhood list
            if Pattern[i] != j:
                Neighborhood.append(Pattern[:i] + j + Pattern[i+1:])
    return Neighborhood

def Neighbors(Pattern, d):
    """
    Finds all neighbors of a given pattern within a specified Hamming distance.

    Parameters:
    Pattern (str): The input pattern string.
    d (int): Maximum Hamming distance.

    Returns:
    list: A list containing all patterns that are within the maximum Hamming distance d of the input pattern.
    """
    # Initialize the neighborhood with the input pattern itself
    Neighborhood = [Pattern]
    # For each distance from 0 to d, find the immediate neighbors of the current patterns
    for j in range(d):
        # For each current neighbor, find its immediate neighbors and add them to the neighborhood
        for Neighbor in Neighborhood:
            Neighborhood = Neighborhood + ImmediateNeighbor(Neighbor)
        # Remove duplicate neighbors to ensure uniqueness in the neighborhood
        Neighborhood = list(set(Neighborhood))
    # Remove any remaining duplicates before returning the final neighborhood
    Neighborhood = list(set(Neighborhood))
    return Neighborhood

def MultiNeighbors(k, d, Dna):
    """
    Calculate the d-neighborhood set of k-mers for each DNA sequence.

    Parameters:
    k -- int, the value of k in k-mers
    d -- int, the maximum distance for the neighborhood
    Dna -- list, a list of input DNA sequences

    Returns:
    Neighborhood -- dict, the d-neighborhood set of k-mers for each DNA sequence
    """
    # Initialize an empty dictionary to store the neighborhood sets for each DNA sequence
    Neighborhood = {}
    # Iterate over each DNA sequence in the input
    for String in Dna:
        # Initialize an empty list for the current DNA sequence's neighborhood set
        Neighborhood[String] = []
        # Iterate over all possible k-mers in the current DNA sequence
        for i in range(len(String) - k + 1):
            # Extract the k-mer at the current position
            Pattern = String[i:i+k]
            # Add the d-neighborhood of the current k-mer to the neighborhood set
            Neighborhood[String] = Neighborhood[String] + Neighbors(Pattern, d)
        # Remove duplicates from the neighborhood set
        Neighborhood[String] = list(set(Neighborhood[String]))
    # Return the neighborhood sets for each DNA sequence
    return Neighborhood

def MotifEnumeration(Dna, k, d):
    """
    Main function: enumerates all common motifs in DNA sequences given certain conditions.

    Parameters:
    Dna - A list containing multiple DNA sequences.
    k - The length of the motif.
    d - The allowed Hamming distance.

    Returns:
    A list containing all motifs that satisfy the conditions.
    """
    # Initialize the motifs list
    Motifs = []
    # Build a dictionary of k-mers and their neighborhoods for each DNA sequence
    Neighborhood = MultiNeighbors(k, d, Dna)
    
    # Iterate through each DNA sequence
    for String in Dna:
        # Iterate through all possible k-mers and their neighborhoods for the current DNA sequence
        for Pattern in Neighborhood[String]:
            # Ensure the current pattern has not been added to the motifs list
            if Pattern not in Motifs:
                # Initialize a counter to count the number of matching DNA sequences
                Count = 0
                # Iterate through each DNA sequence
                for String2 in Dna:
                    # Iterate through all possible k-mers for the current DNA sequence
                    for i in range(len(String2) - k + 1):
                        # If the Hamming distance between the current pattern and the k-mer is no more than d, increment the counter
                        if HammingDist(Pattern, String2[i:i+k]) <= d:
                            Count += 1
                            break
                # If all DNA sequences match, add the current pattern to the motifs list
                if Count == len(Dna):
                    Motifs.append(Pattern)
    
    # Remove duplicate motifs
    Motifs = list(set(Motifs))
    # Return all motifs that satisfy the conditions
    return(Motifs)

# Open the specified text file and prepare to read its data
with open("dataset_30302_8.txt", "r") as f:
    # Read the first line of the file, split it into two parts, and convert them to integers.
    # It is assumed that the first line contains two integers separated by a space, representing k and d.
    k, d = f.readline().split()
    k, d = int(k), int(d)
    # Read the second line of the file and split it into a list of strings.
    # It is assumed that the second line contains DNA sequences separated by spaces.
    Dna = f.readline().split()
    # Call the MotifEnumeration function with the list of DNA sequences, k, and d, to obtain a set of motifs that meet the criteria.
    Motifs = MotifEnumeration(Dna, k, d)
    # Convert the obtained list of motifs into a string separated by spaces and print it out.
    print(k, d)
    print(Dna)
    print(" ".join(Motifs))
    ##############################

import numpy as np

# DNA sequences
sequences = [
    "TCGGGGGTTTTT",
    "CCGGTGACTTAC",
    "ACGGGGATTTTC",
    "TTGGGGACTTTT",
    "AAGGGGACTTCC",
    "TTGGGGACTTCC",
    "TCGGGGATTCAT",
    "TCGGGGATTCCT",
    "TAGGGGAACTAC",
    "TCGGGTATAACC"
]

# Convert sequences into a NumPy array
matrix = np.array([list(seq) for seq in sequences])

# Number of sequences
num_sequences = len(sequences)

# Initialize a dictionary to count nucleotides (A, T, G, C)
counts = {nucleotide: np.sum(matrix == nucleotide, axis=0) for nucleotide in 'ATGC'}

# Normalize the counts to get probabilities
probabilities = {nucleotide: counts[nucleotide] / num_sequences for nucleotide in 'ATGC'}

# Compute entropy for each column
entropy = np.zeros(matrix.shape[1])  # Initialize an array for column entropies

for nucleotide, probs in probabilities.items():
    # Compute entropy only for non-zero probabilities
    entropy += -probs * np.log2(probs, where=(probs > 0))

# Sum entropy across all columns
total_entropy = np.sum(entropy)

# Print the total entropy
print(f"Total Entropy: {total_entropy:.4f}")
################################

import itertools

def HammingDistance(p, q):
    """
    Computes the Hamming distance between two strings of equal length.
    """
    return sum(1 for pi, qi in zip(p, q) if pi != qi)

def DistanceBetweenPatternAndString(Pattern, string):
    """
    Computes the distance between a k-mer Pattern and a string.
    The distance is the minimum Hamming distance between Pattern and any k-mer in the string.
    """
    k = len(Pattern)
    min_distance = float('inf')
    for i in range(len(string) - k + 1):
        substring = string[i:i + k]
        distance = HammingDistance(Pattern, substring)
        if distance < min_distance:
            min_distance = distance
    return min_distance

def DistanceBetweenPatternAndDna(Pattern, Dna):
    """
    Computes the total distance between a k-mer Pattern and a collection of strings Dna.
    """
    return sum(DistanceBetweenPatternAndString(Pattern, string) for string in Dna)

def MedianString(k, Dna):
    """
    Finds a k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers.
    """
    # Generate all possible k-mers
    nucleotides = "ACGT"
    possible_kmers = itertools.product(nucleotides, repeat=k)
    best_pattern = None
    best_distance = float('inf')
    
    for kmer_tuple in possible_kmers:
        kmer = "".join(kmer_tuple)
        distance = DistanceBetweenPatternAndDna(kmer, Dna)
        if distance < best_distance:
            best_distance = distance
            best_pattern = kmer
    
    return best_pattern

# Read the dataset
with open("dataset_30304_9.txt", "r") as file:
    data = file.read().splitlines()
    k = int(data[0])  # First line is the value of k
    Dna = data[1].split()  # Second line contains the DNA strings

# Solve MedianString
median_kmer = MedianString(k, Dna)
print(f"Median k-mer: {median_kmer}")
######################

import numpy as np

# Profile matrix
profile = {
    'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0, 0.0],
    'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6, 0.0],
    'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
    'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.5, 0.8, 0.7, 0.3, 0.4, 0.4],
}

# Sequence
sequence = "TCGTGGATTTCC"

# Calculate probability
probability = 1.0
for i, nucleotide in enumerate(sequence):
    probability *= profile[nucleotide][i]

print(f"Pr({sequence} | Profile) = {probability:.6f}")

#######################################

def ComputeProbability(kmer, profile):
    """
    Compute the probability of a k-mer given a profile matrix.

    Parameters:
        kmer (str): The k-mer string.
        profile (dict): A profile matrix as a dictionary with keys 'A', 'C', 'G', 'T'.

    Returns:
        float: The probability of the k-mer given the profile.
    """
    probability = 1.0
    for i, nucleotide in enumerate(kmer):
        probability *= profile[nucleotide][i]
    return probability

def ProfileMostProbableKmer(text, k, profile):
    """
    Find the Profile-most probable k-mer in a string.

    Parameters:
        text (str): The input text.
        k (int): The length of the k-mer.
        profile (dict): A profile matrix as a dictionary with keys 'A', 'C', 'G', 'T'.

    Returns:
        str: The Profile-most probable k-mer in the text.
    """
    max_probability = -1
    most_probable_kmer = text[:k]  # Initialize with the first k-mer

    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        probability = ComputeProbability(kmer, profile)
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer

    return most_probable_kmer

# Read dataset
with open("dataset_30305_3.txt", "r") as file:
    data = file.read().splitlines()
    text = data[0]  # First line is the text
    k = int(data[1])  # Second line is k
    profile_matrix = [list(map(float, line.split())) for line in data[2:]]  # Profile matrix starts from line 3

# Convert profile_matrix into a dictionary
profile = {
    'A': profile_matrix[0],
    'C': profile_matrix[1],
    'G': profile_matrix[2],
    'T': profile_matrix[3]
}

# Solve the problem
result = ProfileMostProbableKmer(text, k, profile)
print(f"Profile-most probable k-mer: {result}")

######################################################

def HammingDistance(p, q):
    """
    Computes the Hamming distance between two strings of equal length.
    """
    return sum(1 for pi, qi in zip(p, q) if pi != qi)

def DistanceBetweenPatternAndStrings(Pattern, Dna):
    """
    Computes the sum of distances between Pattern and each string in Dna.

    Parameters:
        Pattern (str): The k-mer pattern.
        Dna (list): A list of DNA strings.

    Returns:
        int: The sum of distances between Pattern and each string in Dna.
    """
    k = len(Pattern)
    total_distance = 0

    for Text in Dna:
        min_distance = float('inf')  # Initialize with infinity
        for i in range(len(Text) - k + 1):
            kmer = Text[i:i + k]
            distance = HammingDistance(Pattern, kmer)
            if distance < min_distance:
                min_distance = distance
        total_distance += min_distance

    return total_distance

# Read dataset
def load_dataset(filename):
    """
    Reads the dataset from the specified file.

    Parameters:
        filename (str): The path to the dataset file.

    Returns:
        tuple: A tuple containing the pattern and a list of DNA strings.
    """
    with open(filename, "r") as file:
        data = file.read().splitlines()
        Pattern = data[0]  # First line is the Pattern
        Dna = data[1].split()  # Second line contains space-separated DNA strings
    return Pattern, Dna

# Main execution
if __name__ == "__main__":
    # Replace 'dataset_30312_1.txt' with your actual dataset filename
    dataset_filename = "dataset_30312_1.txt"
    Pattern, Dna = load_dataset(dataset_filename)

    # Solve the problem
    result = DistanceBetweenPatternAndStrings(Pattern, Dna)
    print(f"Distance between Pattern and Dna: {result}")
##################################3
from itertools import product

def HammingDistance(p, q):
    """
    Computes the Hamming distance between two strings of equal length.
    """
    return sum(1 for pi, qi in zip(p, q) if pi != qi)

def Neighbors(Pattern, d):
    """
    Generates all neighbors of a pattern within a given Hamming distance.

    Parameters:
        Pattern (str): The input pattern.
        d (int): Maximum Hamming distance.

    Returns:
        set: A set of all neighbors of the pattern within Hamming distance d.
    """
    if d == 0:
        return {Pattern}
    if len(Pattern) == 1:
        return {'A', 'C', 'G', 'T'}
    
    neighbors = set()
    suffix_neighbors = Neighbors(Pattern[1:], d)
    for text in suffix_neighbors:
        if HammingDistance(Pattern[1:], text) < d:
            for nucleotide in 'ACGT':
                neighbors.add(nucleotide + text)
        else:
            neighbors.add(Pattern[0] + text)
    return neighbors

def MotifEnumeration(Dna, k, d):
    """
    Finds all (k, d)-motifs in a collection of DNA strings.

    Parameters:
        Dna (list): A list of DNA strings.
        k (int): The length of the motif.
        d (int): Maximum Hamming distance.

    Returns:
        list: A list of all (k, d)-motifs.
    """
    Patterns = set()
    for string in Dna:
        for i in range(len(string) - k + 1):
            Pattern = string[i:i+k]
            for neighbor in Neighbors(Pattern, d):
                if all(any(HammingDistance(neighbor, string[j:j+k]) <= d 
                           for j in range(len(string) - k + 1)) for string in Dna):
                    Patterns.add(neighbor)
    return sorted(Patterns)

# Input data
k = 3
d = 1
Dna = ["ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"]

# Solve the problem
result = MotifEnumeration(Dna, k, d)
print(" ".join(result))


########

import numpy as np

# Profile matrix
profile = {
    "A": [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
    "C": [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
    "G": [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
    "T": [0.3, 0.1, 0.0, 0.4, 0.5, 0.0],
}

# Strings to check
strings = [
    "ACGTTA",
    "ACGCGA",
    "ATGCTA",
    "AGGTCA",
    "ACGTTT",
    "AAGTGA",
]

def is_consensus_string(profile, string):
    """
    Check if a string is a consensus string for the given profile matrix.
    """
    for i, nucleotide in enumerate(string):
        # Check if the nucleotide has the maximum probability at its position
        max_prob = max(profile["A"][i], profile["C"][i], profile["G"][i], profile["T"][i])
        if profile[nucleotide][i] != max_prob:
            return False
    return True

# Check each string
consensus_strings = []
for string in strings:
    if is_consensus_string(profile, string):
        consensus_strings.append(string)

# Output the consensus strings
print("Consensus Strings:", consensus_strings)


###############

def hamming_distance(s1, s2):
    """
    Compute the Hamming distance between two strings of equal length.
    """
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def total_distance(candidate, motif_matrix):
    """
    Compute the total Hamming distance between a candidate k-mer and the motif matrix.
    """
    total = 0
    for row in motif_matrix:
        min_distance = float("inf")
        for i in range(len(row) - len(candidate) + 1):
            kmer = row[i:i + len(candidate)]
            min_distance = min(min_distance, hamming_distance(candidate, kmer))
        total += min_distance
    return total

def generate_all_kmers(k):
    """
    Generate all possible k-mers of length k using the DNA alphabet.
    """
    from itertools import product
    return [''.join(kmer) for kmer in product('ACGT', repeat=k)]

def find_median_strings(motif_matrix, k):
    """
    Find the k-mer(s) that minimize the total Hamming distance to the motif matrix.
    """
    all_kmers = generate_all_kmers(k)
    min_distance = float("inf")
    median_strings = []

    for kmer in all_kmers:
        dist = total_distance(kmer, motif_matrix)
        if dist < min_distance:
            min_distance = dist
            median_strings = [kmer]
        elif dist == min_distance:
            median_strings.append(kmer)

    return median_strings

# Input motif matrix
motif_matrix = [
    "CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC",
    "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC",
    "GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG",
]

# Set k (length of k-mer)
k = 7

# Find the median strings
median_strings = find_median_strings(motif_matrix, k)
print("Median Strings:", median_strings)
####################


def compute_probability(kmer, profile):
    """
    Compute the probability of a k-mer given a profile matrix.

    Parameters:
    kmer (str): The k-mer to compute the probability for.
    profile (dict): A dictionary representing the profile matrix. Keys are nucleotides (A, C, G, T),
                    and values are lists of probabilities at each position.

    Returns:
    float: The probability of the k-mer given the profile matrix.
    """
    probability = 1.0
    for i, nucleotide in enumerate(kmer):
        probability *= profile[nucleotide][i]
    return probability

# Profile matrix
profile = {
    "A": [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
    "C": [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
    "G": [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
    "T": [0.3, 0.1, 0.0, 0.4, 0.5, 0.0],
}

# k-mer to compute probability for
kmer = "GAGCTA"

# Compute the probability
probability = compute_probability(kmer, profile)
print(f"Pr({kmer}|Profile) = {probability:.6f}")
#########################

###################################################

















