# Function to compute Hamming distance
def HammingDistance(p, q):
    """
    Computes the Hamming distance between two strings of equal length.
    """
    return sum(1 for pi, qi in zip(p, q) if pi != qi)

# Hamming distance calculation
p = "CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA"
q = "CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG"
hamming_distance = HammingDistance(p, q)
print("Hamming Distance:", hamming_distance)

# Function to compute Skew values
def ComputeSkew(Genome):
    """
    Computes the skew diagram for a given genome string.

    Parameters:
        Genome (str): The DNA string.

    Returns:
        list: A list of integers representing Skew_i values from i = 0 to |Genome|.
    """
    skew = [0]  # Start with Skew_0 = 0
    for nucleotide in Genome:
        if nucleotide == 'G':
            skew.append(skew[-1] + 1)
        elif nucleotide == 'C':
            skew.append(skew[-1] - 1)
        else:
            skew.append(skew[-1])
    return skew

# Skew calculation
Genome = "CATTCCAGTACTTCGATGATGGCGTGAAGA"
skew_values = ComputeSkew(Genome)
min_value = min(skew_values)
min_indices = [i for i, value in enumerate(skew_values) if value == min_value]
print("Minimum Skew Index:", min_indices)

# Function to count pattern occurrences
def PatternCount(Text, Pattern):
    """
    Counts the number of times a Pattern appears in a given Text, including overlapping occurrences.

    Parameters:
        Text (str): The string to search within.
        Pattern (str): The pattern to count.

    Returns:
        int: The number of occurrences of the Pattern in the Text.
    """
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i:i + len(Pattern)] == Pattern:
            count += 1
    return count

# Pattern count
Text = "CGTGACAGTGTATGGGCATCTTT"
Pattern = "TGT"
count = PatternCount(Text, Pattern)
print("Pattern Count:", count)

# Function to compute the d-neighborhood
from itertools import product

def Neighborhood(Pattern, d):
    """
    Finds the d-neighborhood of a k-mer Pattern.

    Parameters:
        Pattern (str): The k-mer to generate neighbors for.
        d (int): The maximum Hamming distance.

    Returns:
        set: A set of k-mers within Hamming distance d of the Pattern.
    """
    nucleotides = ['A', 'C', 'G', 'T']
    neighborhood = set()

    # Generate all possible k-mers of the same length as Pattern
    for kmer in product(nucleotides, repeat=len(Pattern)):
        kmer = ''.join(kmer)
        if HammingDistance(Pattern, kmer) <= d:
            neighborhood.add(kmer)

    return neighborhood

# d-neighborhood calculation
Pattern = "CCAGTCAATG"
d = 1
neighborhood = Neighborhood(Pattern, d)
print("d-Neighborhood Size:", len(neighborhood))
