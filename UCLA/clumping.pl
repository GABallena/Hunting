use strict;
use warnings;
use List::MoreUtils qw(uniq);

# Function to find k-mers forming (L, t)-clumps in Genome
sub ClumpFinding {
    my ($genome, $k, $L, $t) = @_;
    my %clumps;

    # Iterate through all possible starting positions of the window of length L
    for (my $i = 0; $i <= length($genome) - $L; $i++) {
        my %kmer_counts;  # To count occurrences of k-mers in the window

        # For each window, count the frequency of each k-mer (including overlaps)
        for (my $j = $i; $j <= $i + $L - $k; $j++) {
            my $kmer = substr($genome, $j, $k);
            $kmer_counts{$kmer}++;
        }

        # Add k-mers to clumps if their count is greater than or equal to t
        foreach my $kmer (keys %kmer_counts) {
            if ($kmer_counts{$kmer} >= $t) {
                $clumps{$kmer} = 1;
            }
        }
    }

    # Return the distinct clumps as a list of sorted k-mers
    return keys %clumps;
}

# Main function to read the dataset and process the data
sub main {
    # Read the genome from the file
    my $filename = 'E_coli.txt';  # Path to the E. coli genome file
    open my $fh, '<', $filename or die "Could not open file '$filename': $!";
    my $genome = do { local $/; <$fh> };  # Read the entire genome into a single string
    close $fh;

    my $k = 9;  # Length of the k-mer
    my $L = 500;  # Length of the sliding window
    my $t = 3;    # Minimum occurrences for a clump

    # Run the ClumpFinding function on the genome
    my @clumps = ClumpFinding($genome, $k, $L, $t);

    # Print the number of distinct 9-mers forming (500, 3)-clumps
    print "Number of distinct 9-mers forming (500, 3)-clumps: " . scalar(@clumps) . "\n";
}

# Run the main function
main();
