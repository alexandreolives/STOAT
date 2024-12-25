import argparse

def extract_vcf_header(input_vcf):
    header_lines = []
    column_header = None
    num_samples = 0

    # Open the VCF file and process line by line
    with open(input_vcf, 'r') as vcf:
        for line in vcf:
            if line.startswith("##"):
                # Replace specific headers with custom definitions only once
                if line.startswith("##INFO=<ID=P,"):
                    continue
                elif line.startswith("##FORMAT=<ID=GT,"):
                    continue
                elif line.startswith("##FILTER="):
                    continue
                header_lines.append(line)  # Collect unmodified metadata lines

            elif line.startswith("#CHROM"):
                column_header = line  # Identify the column header line
                # Add new or updated headers before the column header
                header_lines.extend([
                    '##INFO=<ID=P,Number=1,Type=Float,Description="P-value from GWAS">',
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                    '##FILTER=<ID=LOWQ,Description="Low quality based on p-value threshold">',
                ])
                # Determine the number of samples
                num_samples = len(line.strip().split('\t')) - 9  # Columns after the FORMAT field
                break  # No need to read further; the header ends here

    return header_lines, column_header, num_samples

def create_vcf_from_gwas(gwas_file, input_vcf, output_vcf):
    # Extract the header lines, column header, and number of samples
    header_lines, column_header, num_samples = extract_vcf_header(input_vcf)

    # Open the output VCF for writing
    with open(output_vcf, 'w') as out_vcf:
        out_vcf.writelines(header_lines)
        out_vcf.write(column_header)

        # Generate placeholder sample data (e.g., "." repeated for each sample)
        sample_placeholder = "\t".join(["."] * num_samples)

        # Process GWAS file and generate VCF body
        with open(gwas_file, 'r') as gwas:
            next(gwas)
            for line in gwas:
                fields = line.strip().split('\t')
                chrom, pos, snarl_id, _, ref, alt, p_value = fields

                # Create placeholder fields
                qual = "."
                filter_field = "PASS" if float(p_value) <= 0.05 else "LOWQ"
                info_field = f"P={p_value}"
                format_field = "GT"

                # Create and write the VCF line
                vcf_line = f"{chrom}\t{pos}\t{snarl_id}\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info_field}\t{format_field}\t{sample_placeholder}\n"
                out_vcf.write(vcf_line)

def main():
    parser = argparse.ArgumentParser(description="Convert GWAS output to a VCF file using an input VCF header.")
    parser.add_argument('-g', '--gwas', type=str, required=True, help="Path to the GWAS output file (TSV format).")
    parser.add_argument('-i', '--input_vcf', type=str, required=True, help="Path to the input VCF file (to extract header information).")
    parser.add_argument('-o', '--output_vcf', type=str, required=True, help="Path to the output VCF file.")
    
    args = parser.parse_args()

    # Run the function with the parsed arguments
    create_vcf_from_gwas(args.gwas, args.input_vcf, args.output_vcf)

if __name__ == "__main__":
    main()

# (base) mbagarre@irsd0440:~/Bureau/STOAT$ bcftools norm tests/vcf_test.vcf -m -any -f ../droso_data/fly/ref.fa > vcf_test.normalized.vcf
# [E::vcf_hdr_read] No sample line
# Failed to read from tests/vcf_test.vcf: could not parse header