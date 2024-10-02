import argparse

# Function to parse the VCF file and create a dictionary with ID as key and POS as value
def parse_vcf_to_dict(vcf_file):
    vcf_dict = {}

    with open(vcf_file, 'r') as f:
        for line in f:
            # Skip header lines
            if line.startswith('#'):
                continue

            # Split the line by tabs
            fields = line.strip().split('\t')
            
            # Extract the POS and ID
            chr = fields[0]  # 2nd column is POS
            pos = fields[1]  # 2nd column is POS
            snarl_id = fields[2]  # 3rd column is ID
            
            # Add to dictionary
            vcf_dict[snarl_id] = (pos, chr)

    return vcf_dict

# Function to match the snarl list against the VCF dictionary
def match_snarl_to_vcf(snarl_file, vcf_dict):
    with open(snarl_file, 'r') as f:
        for line in f:
            # Parse each snarl_id and modify it (replace _ with >)
            original_snarl_id = line.strip().split('\t')[0]
            modified_snarl_id = f">{original_snarl_id.replace('_', '>')}"

            # Also create the reversed version of the snarl ID (for cases like >25>12)
            reversed_snarl_id = f">{original_snarl_id.split('_')[1]}>{original_snarl_id.split('_')[0]}"

            # Check if the modified snarl ID or its reversed form exists in the VCF dictionary
            if modified_snarl_id in vcf_dict:
                print(f"{vcf_dict[modified_snarl_id][1]} {vcf_dict[modified_snarl_id][0]}")
            elif reversed_snarl_id in vcf_dict:
                print(f"{vcf_dict[reversed_snarl_id][1]} {vcf_dict[reversed_snarl_id][0]}")
            else:
                print(f"No match found for {modified_snarl_id} or {reversed_snarl_id}")

def main():
    parser = argparse.ArgumentParser(description="Parse VCF and snarl files, then match snarl IDs with positions.")
    parser.add_argument('--vcf', required=True, help="Path to the VCF file")
    parser.add_argument('--snarl', required=True, help="Path to the snarl ID file")
    args = parser.parse_args()
  
    vcf_dict = parse_vcf_to_dict(args.vcf)
    match_snarl_to_vcf(args.snarl, vcf_dict)

"""
  python3 src/snarl_annotation.py --vcf ../../droso_data/vcf/merged_output.vcf --snarl output/droso_significative_snarl.tsv > output/annotation_droso.txt
"""

if __name__ == "__main__":
    main()
