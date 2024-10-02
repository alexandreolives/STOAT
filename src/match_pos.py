import argparse

# Function to read positions from a file
def read_positions(file_name):
    positions = []
    with open(file_name, 'r') as file:
        for line in file:
            chrom, pos = line.strip().split()
            positions.append((chrom, int(pos)))  # Store as tuple (chromosome, position)
    return positions

# Function to find matches within ±1000 positions
def find_matches(positions1, positions2):
    matches = []

    for chrom1, pos1 in positions1:
        for chrom2, pos2 in positions2:
            if chrom1 == chrom2 and abs(pos1 - pos2) <= 1000:
                matches.append((chrom1, pos1, pos2))
                break
    return matches

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Match positions from two files within ±1000 positions.")
    parser.add_argument('file1', type=str, help='The first input file containing positions.')
    parser.add_argument('file2', type=str, help='The second input file containing positions.')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Read positions from both files
    positions1 = read_positions(args.file1)
    positions2 = read_positions(args.file2)

    # Find matches
    matches = find_matches(positions1, positions2)

    # Output matches
    for chrom, pos1, pos2 in matches:
        print(f"Chromosome: {chrom}, Position 1: {pos1}, Position 2: {pos2}")

if __name__ == '__main__':
    main()
