def PatternMatching(Pattern, Genome):
    """
    Finds all starting positions where Pattern appears as a substring of Genome.

    Parameters:
        Pattern (str): The substring to search for.
        Genome (str): The string in which to search for Pattern.

    Returns:
        str: A space-separated string of all starting positions (0-based indexing).
    """
    positions = []
    pattern_length = len(Pattern)
    
    # Find all starting positions of Pattern in Genome
    for i in range(len(Genome) - pattern_length + 1):
        if Genome[i:i + pattern_length] == Pattern:
            positions.append(i)
    
    # Return positions as a space-separated string
    return ' '.join(map(str, positions))

# Main function
def main():
    # File containing the Vibrio cholerae genome
    genome_filename = 'Vibrio_cholerae.txt'

    # Pattern to search for
    Pattern = "CTTGATCAT"

    # Read the genome file
    with open(genome_filename, 'r') as file:
        Genome = file.read().strip()

    # Find all starting positions of the pattern
    result = PatternMatching(Pattern, Genome)

    # Print the result
    print(result)

# Run the main function
if __name__ == "__main__":
    main()
