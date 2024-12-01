from collections import Counter
import random

# Step 1: Parse FASTA file
def parse_fasta(file_path):
    """Parses a FASTA file and returns a list of sequences."""
    with open(file_path, 'r') as file:
        sequences = []
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

# Step 2: GibbsSampler
def profile_with_pseudocounts(motifs, k):
    t = len(motifs)
    profile = {base: [1] * k for base in 'ACGT'}
    for motif in motifs:
        for i, base in enumerate(motif):
            profile[base][i] += 1
    for base in profile:
        profile[base] = [count / (t + 4) for count in profile[base]]
    return profile

def profile_random_kmer(dna, k, profile):
    n = len(dna)
    kmers = [dna[i:i + k] for i in range(n - k + 1)]
    probs = []
    for kmer in kmers:
        prob = 1
        for i, base in enumerate(kmer):
            prob *= profile[base][i]
        probs.append(prob)
    total_prob = sum(probs)
    probs = [p / total_prob for p in probs]
    return random.choices(kmers, weights=probs, k=1)[0]

def score(motifs):
    """Calculates the score of a set of motifs."""
    consensus = ''.join(Counter(col).most_common(1)[0][0] for col in zip(*motifs))
    return sum(len(motifs) - Counter(col)[consensus[i]] for i, col in enumerate(zip(*motifs)))

def gibbs_sampler(dna_list, k, t, N):
    """Runs the Gibbs Sampler algorithm."""
    motifs = [dna[random.randint(0, len(dna) - k):][:k] for dna in dna_list]
    best_motifs = motifs[:]
    best_score = score(best_motifs)

    for _ in range(N):
        i = random.randint(0, t - 1)
        motifs_minus_one = [motifs[j] for j in range(len(motifs)) if j != i]
        profile = profile_with_pseudocounts(motifs_minus_one, k)
        motifs[i] = profile_random_kmer(dna_list[i], k, profile)
        current_score = score(motifs)
        if current_score < best_score:
            best_motifs = motifs[:]
            best_score = current_score

    return best_motifs

def gibbs_sampler_with_random_starts(dna_list, k, t, N, random_starts=20):
    """Runs the Gibbs Sampler algorithm with multiple random starts."""
    best_motifs = None
    best_score = float('inf')

    for _ in range(random_starts):
        motifs = gibbs_sampler(dna_list, k, t, N)
        current_score = score(motifs)
        if current_score < best_score:
            best_motifs = motifs[:]
            best_score = current_score

    return best_motifs

# Main function
if __name__ == "__main__":
    # File name is adjusted to upstream250.fasta
    file_path = "upstream250.fasta"
    
    # Step 1: Read sequences from the FASTA file
    sequences = parse_fasta(file_path)
    
    # Step 2: Set parameters for GibbsSampler
    k = 20  # Motif length
    t = len(sequences)  # Number of sequences
    N = 100  # Number of iterations for Gibbs sampling

    # Step 3: Run GibbsSampler
    best_motifs = gibbs_sampler_with_random_starts(sequences, k, t, N, random_starts=20)

    # Step 4: Output the results
    print("\nBest Motifs:")
    print("\n".join(best_motifs))
