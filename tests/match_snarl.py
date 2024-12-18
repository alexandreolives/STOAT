import argparse
import re

def extract_integers(s):
    # Find all sequences of digits in the string
    return [int(num) for num in re.findall(r'\d+', s)]

def identify_header(f):
    """
    Identify indices for CHR, POS, and P in the header line.
    Assumes the header contains these columns.
    """
    header = f.readline().strip().split()
    return header.index("SNARL"), header.index("P")

# Define the file parsing function
def parse_file(file_path, threshold=1.0):
    """
    Parses a file containing pairs of numbers separated by an underscore.
    Returns a set of tuples for easier comparison.
    """
    parsed = set()
    total_count = 0
    not_count = 0
    with open(file_path, 'r') as f:
        snarl_idx, p_idx = identify_header(f)
        for line in f : 
            fields = line.strip().split()
            snarl = fields[snarl_idx]
            if "_" in snarl : 
                snarl_tuple = tuple(snarl.split('_'))
            else :
                snarl_tuple = tuple(extract_integers(snarl))
            p = fields[p_idx]
            if p > threshold :
                not_count += 1
                continue
            parsed.add(snarl_tuple)
            total_count += 1  

    return parsed, total_count, not_count


def compare(file1_pairs, file2_pairs):
    not_found = []
    match_count = 0
    not_match_count = 0
    for snarl in file1_pairs:
        # Convert pair to tuple to ensure order-insensitive comparison
        direct = f'{snarl[0]}_{snarl[1]}'
        reversed = f'{snarl[11]}_{snarl[0]}' 
        if direct not in file2_pairs and reversed not in file2_pairs:
            not_found.append(direct)  # Join the pair back to string format
            not_match_count += 1
            print(not_found)
        else :
            match_count += 1
    return not_found, not_match_count, match_count

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Compare pairs from two files and output not found pairs.")
    parser.add_argument(
        "-f1", "--file1", required=True, help="Path to the first input file (e.g., 'splited_snarl_plink.txt')."
    )
    parser.add_argument(
        "-f2", "--file2", required=True, help="Path to the second input file (e.g., 'snarl.snarl.sorted.tsv')."
    )
    parser.add_argument('-s', '--significativity', type=float, required=False, help='P-value threshold for significance filtering')
    parser.add_argument(
        "-o", "--output", required=False, help="Path to the output file (e.g., 'snarl_not_found.txt')."
    )
    args = parser.parse_args()
    significativity = args.significativity if args.significativity else 1.0

    output = args.output if args.output else "snarl_not_found.tsv"
    
    # Parse the input files
    dict_file1, ref_count, ref_not_count= parse_file(args.file1, significativity)
    dict_file2, input_count, input_not_count = parse_file(args.file2, significativity)
    
    not_found, not_match_count, match_count = compare(dict_file1, dict_file2)

    # Compare file2 with the dictionary and count matches with the given threshold
    print(f"Number of unique variant significatif in {args.file1} : ", ref_count, "/", ref_not_count)
    print(f"Number of unique variant significatif in {args.file2} : ", input_count,  "/", input_not_count)
    print(f"Number of variant not_match_count : ", not_match_count)
    print(f"Number of matching pairs with p-value < {significativity}: {match_count}")

    # Write the not found pairs to the output file
    with open(output, 'w') as output_file:
        output_file.write('\n'.join(not_found))
    
    print(f"Processing complete. Not found pairs written to '{output}'.")

if __name__ == "__main__":
    main()

# INPUT : one column with only snarl per file 

# python3 match_snarl.py -f1 plink.normalized.assoc.tsv -f2 snarl.normalized.tsv 