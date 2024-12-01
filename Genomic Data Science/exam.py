from Bio import SeqIO

# Helper functions
def find_orfs(sequence, frame):
    """Find all ORFs in a given reading frame."""
    orfs = []
    sequence = sequence[frame-1:]  # Adjust frame (1, 2, or 3) to 0-based index
    start = None
    for i in range(len(sequence) - 2):
        codon = sequence[i:i+3]
        if codon == 'ATG' and start is None:
            start = i
        if codon in ['TAA', 'TAG', 'TGA'] and start is not None:
            orfs.append(sequence[start:i+3])
            start = None
    return orfs


def parse_fasta(file_name):
    """Parse the FASTA file and return sequences and headers."""
    sequences = []
    headers = []
    for record in SeqIO.parse(file_name, "fasta"):
        headers.append(record.id)
        sequences.append(str(record.seq))
    return headers, sequences


def get_repeat_counts(sequences, n):
    """Get the count of repeats of length n."""
    repeat_counts = {}
    for seq in sequences:
        for i in range(len(seq) - n + 1):
            repeat = seq[i:i+n]
            if repeat in repeat_counts:
                repeat_counts[repeat] += 1
            else:
                repeat_counts[repeat] = 1
    return repeat_counts


# 1. Number of records in the FASTA file
def num_records(file_name):
    headers, _ = parse_fasta(file_name)
    return len(headers)


# 2. Length of the longest sequence
def longest_sequence_length(file_name):
    _, sequences = parse_fasta(file_name)
    return max(len(seq) for seq in sequences)


# 3. Length of the shortest sequence
def shortest_sequence_length(file_name):
    _, sequences = parse_fasta(file_name)
    return min(len(seq) for seq in sequences)


# 4. Length of the longest ORF in reading frame 2
def longest_orf_frame_2(file_name):
    _, sequences = parse_fasta(file_name)
    longest_orf = 0
    for seq in sequences:
        orfs = find_orfs(seq, 2)
        for orf in orfs:
            longest_orf = max(longest_orf, len(orf))
    return longest_orf


# 5. Starting position of the longest ORF in reading frame 3
def longest_orf_start_frame_3(file_name):
    _, sequences = parse_fasta(file_name)
    longest_orf = 0
    start_pos = 0
    for seq in sequences:
        orfs = find_orfs(seq, 3)
        for orf in orfs:
            if len(orf) > longest_orf:
                longest_orf = len(orf)
                start_pos = seq.find(orf) + 1  # 1-based indexing
    return start_pos


# 6. Length of the longest ORF in any forward reading frame
def longest_orf_any_frame(file_name):
    _, sequences = parse_fasta(file_name)
    longest_orf = 0
    for seq in sequences:
        for frame in range(1, 4):
            orfs = find_orfs(seq, frame)
            for orf in orfs:
                longest_orf = max(longest_orf, len(orf))
    return longest_orf


# 7. Length of the longest forward ORF in the sequence with the identifier gi|142022655|gb|EQ086233.1|16
def longest_orf_in_specific_seq(file_name, identifier):
    headers, sequences = parse_fasta(file_name)
    if identifier in headers:
        idx = headers.index(identifier)
        longest_orf = 0
        for frame in range(1, 4):
            orfs = find_orfs(sequences[idx], frame)
            for orf in orfs:
                longest_orf = max(longest_orf, len(orf))
        return longest_orf
    else:
        return None


# 8. Most frequent repeat of length 6
def most_frequent_repeat(file_name, n=6):
    _, sequences = parse_fasta(file_name)
    repeat_counts = get_repeat_counts(sequences, n)
    most_frequent = max(repeat_counts, key=repeat_counts.get)
    return repeat_counts[most_frequent]


# 9. Number of different repeats of length 12 occurring Max times
def max_repeats_length_12(file_name, n=12):
    _, sequences = parse_fasta(file_name)
    repeat_counts = get_repeat_counts(sequences, n)
    max_occurrences = max(repeat_counts.values())
    return len([repeat for repeat in repeat_counts if repeat_counts[repeat] == max_occurrences])


# 10. Most frequent repeat of length 7
def most_frequent_repeat_length_7(file_name, n=7):
    _, sequences = parse_fasta(file_name)
    repeat_counts = get_repeat_counts(sequences, n)
    most_frequent = max(repeat_counts, key=repeat_counts.get)
    return most_frequent


# Example usage
file_name = "dna2.fasta"

# Answers for the exam
print(f"1. Number of records: {num_records(file_name)}")
print(f"2. Length of the longest sequence: {longest_sequence_length(file_name)}")
print(f"3. Length of the shortest sequence: {shortest_sequence_length(file_name)}")
print(f"4. Length of the longest ORF in reading frame 2: {longest_orf_frame_2(file_name)}")
print(f"5. Starting position of the longest ORF in reading frame 3: {longest_orf_start_frame_3(file_name)}")
print(f"6. Length of the longest ORF in any forward reading frame: {longest_orf_any_frame(file_name)}")
print(f"7. Length of the longest forward ORF in the sequence with identifier gi|142022655|gb|EQ086233.1|16: {longest_orf_in_specific_seq(file_name, 'gi|142022655|gb|EQ086233.1|16')}")
print(f"8. Most frequent repeat of length 6: {most_frequent_repeat(file_name, 6)}")
print(f"9. Number of repeats of length 12 occurring Max times: {max_repeats_length_12(file_name, 12)}")
print(f"10. Most frequent repeat of length 7: {most_frequent_repeat_length_7(file_name, 7)}")
#######################3

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Function to read and parse the lambda virus genome from a FASTA file
def read_genome(filename):
    """Reads the lambda virus genome from a FASTA file."""
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # Ignore header lines (those starting with '>')
            if not line.startswith('>'):
                genome += line.rstrip()  # Concatenate all sequence lines
    return genome

# Function to find occurrences of a pattern (exact matches) in the text (genome)
def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i + j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

# Function to find the reverse complement of a DNA sequence
def reverse_complement(s):
    """Returns the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(s))

# Function to read a FASTQ file
def read_fastq(filename):
    """Read and parse a FASTQ file."""
    reads = []
    with open(filename, 'r') as file:
        for record in SeqIO.parse(file, "fastq"):
            reads.append(record)
    return reads

# Function to find poor-quality cycle in the FASTQ reads
def find_poor_quality_cycle(fastq_reads):
    """Find the sequencing cycle with the poorest average quality."""
    quality_scores = np.array([list(record.letter_annotations["phred_quality"]) for record in fastq_reads]).T
    average_quality = np.mean(quality_scores, axis=1)
    poor_quality_cycle = np.argmin(average_quality)  # Find the index of the minimum average quality
    return poor_quality_cycle, average_quality[poor_quality_cycle]

# Function for approximate matching (with up to 2 mismatches)
def naive_2mm(p, t, max_mismatches=2):
    """Find approximate matches of a pattern in a text, allowing up to 2 mismatches."""
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        mismatches = 0
        for j in range(len(p)):  # loop over characters
            if t[i + j] != p[j]:  # compare characters
                mismatches += 1
            if mismatches > max_mismatches:
                break
        if mismatches <= max_mismatches:
            occurrences.append(i)  # Record position if mismatches are within limit
    return occurrences

# Read the lambda virus genome
genome = read_genome('lambda_virus.fa')

# Example Pattern Matching - AGGT and its reverse complement (ACCT)
pattern1 = "AGGT"
reverse_pattern1 = reverse_complement(pattern1)
occurrences_aggt = naive(pattern1, genome)
occurrences_acct = naive(reverse_pattern1, genome)

# Total occurrences of AGGT and its reverse complement
total_occurrences = len(occurrences_aggt) + len(occurrences_acct)
print(f"Total occurrences of AGGT and its reverse complement (ACCT): {total_occurrences}")

# Example Pattern Matching - TTAA and its reverse complement (same in this case)
pattern2 = "TTAA"
reverse_pattern2 = reverse_complement(pattern2)
occurrences_ttaa = naive(pattern2, genome)

# Number of occurrences of TTAA and its reverse complement (same in this case)
print(f"Occurrences of TTAA (reverse complement is the same): {len(occurrences_ttaa)}")

# Leftmost occurrence of a pattern ACTAAGT or its reverse complement
pattern3 = "ACTAAGT"
reverse_pattern3 = reverse_complement(pattern3)
leftmost_occurrence_actaagt = min(naive(pattern3, genome) + naive(reverse_pattern3, genome))
print(f"Leftmost occurrence of ACTAAGT or its reverse complement is at offset: {leftmost_occurrence_actaagt}")

# Leftmost occurrence of AGTCGA or its reverse complement
pattern4 = "AGTCGA"
reverse_pattern4 = reverse_complement(pattern4)
leftmost_occurrence_agtcga = min(naive(pattern4, genome) + naive(reverse_pattern4, genome))
print(f"Leftmost occurrence of AGTCGA or its reverse complement is at offset: {leftmost_occurrence_agtcga}")

# Now let's analyze the FASTQ file with real DNA sequencing reads
fastq_reads = read_fastq("ERR037900_1.first1000.fastq")
print(f"Total number of reads in FASTQ file: {len(fastq_reads)}")

# Find the poor quality cycle in the sequencing reads
poor_cycle, poor_quality = find_poor_quality_cycle(fastq_reads)
print(f"The poor quality cycle is at position {poor_cycle} with an average quality score of {poor_quality:.2f}.")

# Approximate matching with up to 2 mismatches for TTCAAGCC in the lambda virus genome
pattern5 = "TTCAAGCC"
approx_occurrences = naive_2mm(pattern5, genome)
print(f"Occurrences of {pattern5} (allowing up to 2 mismatches): {len(approx_occurrences)}")

# Leftmost occurrence of AGGAGGTT with up to 2 mismatches in the lambda virus genome
pattern6 = "AGGAGGTT"
approx_leftmost_occurrence = min(naive_2mm(pattern6, genome))
print(f"Leftmost occurrence of {pattern6} (allowing up to 2 mismatches) is at offset: {approx_leftmost_occurrence}")
##################


import bisect

# ----------------------------------------------
# Naive Exact Matching Algorithm
# ----------------------------------------------
def naive_exact_matching(text, pattern):
    n = len(text)
    m = len(pattern)
    alignments = 0
    comparisons = 0
    
    for i in range(n - m + 1):  # for each possible alignment
        alignments += 1
        for j in range(m):  # compare the pattern to the substring
            comparisons += 1
            if text[i + j] != pattern[j]:
                break
    return alignments, comparisons

# Example usage:
text = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
pattern = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"

alignments, comparisons = naive_exact_matching(text, pattern)
print(f"Naive Exact Matches: {alignments} alignments, {comparisons} comparisons")

# ----------------------------------------------
# Boyer-Moore Algorithm
# ----------------------------------------------
def build_bad_character_table(pattern):
    m = len(pattern)
    table = {}
    for i in range(m - 1):
        table[pattern[i]] = m - i - 1
    return table

def boyer_moore(text, pattern):
    n = len(text)
    m = len(pattern)
    table = build_bad_character_table(pattern)
    alignments = 0
    comparisons = 0
    i = m - 1  # Start with the last character of the pattern

    while i < n:
        j = m - 1  # Pattern starts from the last character
        while j >= 0 and text[i] == pattern[j]:
            comparisons += 1
            i -= 1
            j -= 1
        
        if j < 0:  # A match is found
            alignments += 1
            i += m  # Shift pattern to the right
        else:
            comparisons += 1
            bad_char_shift = table.get(text[i], m)  # Default shift if no bad character is found
            i += max(bad_char_shift, m - j - 1)  # Shift pattern based on the bad character rule
    
    return alignments, comparisons

# Example usage:
alignments, comparisons = boyer_moore(text, pattern)
print(f"Boyer-Moore Matches: {alignments} alignments, {comparisons} comparisons")

# ----------------------------------------------
# K-mer Index for Approximate Matching
# ----------------------------------------------
class Index:
    def __init__(self, t, k):
        self.k = k  # k-mer length
        self.index = []
        for i in range(len(t) - k + 1):  # For each k-mer
            self.index.append((t[i:i + k], i))
        self.index.sort()  # Sort k-mers

    def query(self, p):
        kmer = p[:self.k]  # Query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # Binary search
        hits = []
        while i < len(self.index):
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def approximate_matching_with_index(t, p, k=8, max_mismatches=2):
    index = Index(t, k)
    matches = set()  # To avoid duplicate matches
    n = len(t)
    m = len(p)
    
    # Sliding window over the text
    for i in range(n - m + 1):
        for j in range(m - k + 1):  # Checking for k-mer matches
            hits = index.query(p[j:j + k])
            for hit in hits:
                # Check if the match has 2 or fewer mismatches
                if sum(1 for a, b in zip(t[hit:hit + m], p) if a != b) <= max_mismatches:
                    matches.add(hit)  # Store match position
    return len(matches)

# Example usage:
matches = approximate_matching_with_index(text, pattern)
print(f"Approximate Matches (with 2 mismatches): {matches} matches")

# ----------------------------------------------
# Total Index Hits for Approximate Matching
# ----------------------------------------------
def total_index_hits(t, p, k=8):
    index = Index(t, k)
    hits = 0
    
    for i in range(len(p) - k + 1):
        hits += len(index.query(p[i:i + k]))
    
    return hits

# Example usage:
hits = total_index_hits(text, pattern)
print(f"Total Index Hits: {hits} hits")

# ----------------------------------------------
# SubseqIndex for Approximate Matching
# ----------------------------------------------
class SubseqIndex:
    def __init__(self, t, k, ival):
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # interval for subsequences
        self.index = []
        for i in range(0, len(t) - k + 1, ival):  # Generate subsequences
            self.index.append((t[i:i + k], i))
        self.index.sort()

    def query(self, p):
        kmer = p[:self.k]  # Query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # Binary search
        hits = []
        while i < len(self.index):
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def approximate_matching_with_subseq_index(t, p, k=8, ival=3, max_mismatches=2):
    index = SubseqIndex(t, k, ival)
    matches = set()  # To avoid duplicate matches
    n = len(t)
    m = len(p)
    
    # Sliding window over the text
    for i in range(n - m + 1):
        for j in range(m - k + 1):  # Checking for subsequence matches
            hits = index.query(p[j:j + k])
            for hit in hits:
                # Check if the match has 2 or fewer mismatches
                if sum(1 for a, b in zip(t[hit:hit + m], p) if a != b) <= max_mismatches:
                    matches.add(hit)  # Store match position
    return len(matches)

# Example usage:
subseq_matches = approximate_matching_with_subseq_index(text, pattern)
print(f"Subseq Approximate Matches (with 2 mismatches): {subseq_matches} matches")
