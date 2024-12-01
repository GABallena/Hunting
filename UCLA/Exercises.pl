use strict;
use warnings;

# Function to find k-mers forming (L, t)-clumps in Genome
sub ClumpFinding {
    my ($genome, $k, $L, $t) = @_;
    my %clumps;

    # Count k-mers in the initial window
    my %kmer_counts;
    for (my $i = 0; $i <= $L - $k; $i++) {
        my $kmer = substr($genome, $i, $k);
        $kmer_counts{$kmer}++;
    }

    # Store the k-mers that form clumps in the initial window
    foreach my $kmer (keys %kmer_counts) {
        $clumps{$kmer} = 1 if $kmer_counts{$kmer} >= $t;
    }

    # Slide the window through the genome
    for (my $i = 1; $i <= length($genome) - $L; $i++) {
        # Remove the k-mer that is sliding out of the window
        my $out_kmer = substr($genome, $i - 1, $k);
        $kmer_counts{$out_kmer}--;

        # Add the new k-mer that is sliding into the window
        my $in_kmer = substr($genome, $i + $L - $k, $k);
        $kmer_counts{$in_kmer}++;

        # Add k-mers to clumps if their count is greater than or equal to t
        foreach my $kmer (keys %kmer_counts) {
            $clumps{$kmer} = 1 if $kmer_counts{$kmer} >= $t;
        }
    }

    # Return the distinct k-mers forming clumps
    return keys %clumps;
}

# Main function to read the dataset and process the data
sub main {
    # Read the genome from the file
    my $filename = 'E_coli.txt';  # Path to the E. coli genome file
    open my $fh, '<', $filename or die "Could not open file '$filename': $!";
    my $genome = do { local $/; <$fh> };  # Read the entire genome as a single string
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
