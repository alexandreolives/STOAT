import argparse
import os 

def identify_header(f):
    """
    Identify indices for CHR, POS, and P in the header line.
    Assumes the header contains these columns.
    """
    header = f.readline().strip().split()
    return header.index("CHR"), header.index("POS"), header.index("P")

# Function to load data from file with significance filtering
def load_file(file, threshold=None) :
    file_dict = {}
    total_succes = 0  # Track total significant entries
    total_fail = 0
    with open(file, 'r') as f:
        #chrom_idx, pos_idx, pvalue_idx = identify_header(f)
        chrom_idx, pos_idx, pvalue_idx = identify_header(f)

        for line in f:
            fields = line.strip().split()
            chrom = fields[chrom_idx]
            pos = fields[pos_idx]
            pvalue = fields[pvalue_idx]

            if threshold :
                if str(pvalue) in ["NA", "nan"] or float(pvalue) > float(threshold) :
                    total_fail += 1
                    continue

            if pos :
                if chrom not in file_dict:
                    file_dict[chrom] = set()
                # Only increment total_succes if pos is new in file1_dict[chrom]
                if pos not in file_dict[chrom]:
                    if "," in pos : # Case mult position 
                        list_pos = pos.split(',')
                        for small_pos in list_pos :
                            file_dict[chrom].add(small_pos)
                            total_succes += 1  
                    else :
                        file_dict[chrom].add(pos)
                        total_succes += 1

    return file_dict, total_succes, total_fail

def compare_dict(file_dict_1, file_dict_2, output_file):
    match_count = 0
    total_count = 0 
    not_match_count = 0

    with open(output_file, 'w') as out:
        out.write("CHR\tPOS\n")  # Write header to output file

        for chrom, positions in file_dict_1.items():
            if chrom in file_dict_2:
                for pos in positions:
                    total_count += 1
                    if pos in file_dict_2[chrom]:
                        out.write(f"{chrom}\t{pos}\n")
                        match_count += 1
                    else:
                        print("chrom/pos : ", chrom, pos)
                        not_match_count += 1
            else: # case not chr find 
                not_match_count += len(positions) 
                total_count += len(positions)

    return not_match_count, match_count, total_count

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Match chromosome and position pairs between two files with significance filtering.')
    
    # Arguments for the input files
    parser.add_argument('ref', type=str, help='Path to the truth input file')
    parser.add_argument('input', type=str, help='Path to the second input file')
    parser.add_argument('-s', '--significativity', type=float, required=False, help='P-value threshold for significance filtering')
    parser.add_argument('-o', '--output', type=str, required=False, help='Output file containing variants that match between files (with p-values < threshold)')

    # Parse the arguments
    args = parser.parse_args()

    significativity = args.significativity if args.significativity else None

    if args.output:
        output = args.output
    else:
        truth_name = os.path.basename(args.ref).replace('.tsv', '')
        input_name = os.path.basename(args.input).replace('.tsv', '')
        output = f"matching_position_{truth_name}_{input_name}_p{significativity}.tsv"

    # Load file1 into a dictionary with the given threshold
    file1_dict, ref_count, ref_fail = load_file(args.ref, significativity)
    file2_dict, input_count, input_fail = load_file(args.input, significativity)

    # Compare file2 with the dictionary and count matches with the given threshold
    not_match_count, match_count, total_count = compare_dict(file1_dict, file2_dict, output)
    print(f"Number of unique variant significatif in {args.ref} : ", ref_count, "not significatif or NA : ", ref_fail)
    print(f"Number of unique variant significatif in {args.input} : ", input_count, "not significatif or NA : ", input_fail)
    print(f"Number of variant not_match_count : ", not_match_count)

    # Calculate and print percentage of matches
    if total_count > 0:
        match_percentage = (match_count / ref_count) * 100

        print(f"Number of matching chromosome and position pairs with p-value < {significativity}: {match_count}")
        print(f"Percentage of matching chr and pos on {truth_name} by {input_name} with p-value < {significativity}: {match_percentage:.2f}%")
    else:
        print("No entries in file2 met the significance threshold.")

    print(f"Output file created at: {os.path.abspath(output)}")

# python3 match_position.py <truth/reference path file> <input path file> -s <threshold>