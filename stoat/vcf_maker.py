import argparse

def extract_vcf_header(input_vcf):
    header_lines = [
        '##fileformat=VCFv4.2\n',
        '##INFO=<ID=P,Number=1,Type=String,Description="P-value">\n',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    ]
    column_header = None
    sample_names = []
    chr_name = []

    # Open the VCF file and process line by line
    with open(input_vcf, 'r') as vcf:
        for line in vcf:
            if line.startswith("##contig=<ID="):
                chr_name.append(line.strip())
            
            elif line.startswith("#CHROM"):
                column_header = line.strip()  # Get the column header line
                columns = column_header.split('\t')
                sample_names = columns[9:]  # Extract sample names from column 10 onward
                break  # No need to read further; the header ends here

    # Add sample names to the header lines
    header_lines.extend([f"{chr}\n" for chr in chr_name])
    header_lines.extend([f"##SAMPLE=<ID={sample}>\n" for sample in sample_names])

    return header_lines, column_header, len(sample_names)

def create_vcf_from_gwas(gwas_file, input_vcf, output_vcf):
    # Extract the header lines, column header, and number of samples
    header_lines, column_header, num_samples = extract_vcf_header(input_vcf)

    # Open the output VCF for writing
    with open(output_vcf, 'w') as out_vcf:
        out_vcf.writelines(header_lines)
        out_vcf.write(column_header + '\n')

        # Process GWAS file and generate VCF body
        with open(gwas_file, 'r') as gwas:
            next(gwas)
            for line in gwas:
                fields = line.strip().split('\t')
                #CHR POS SNARL TYPE	REF	ALT	P_FISHER
                chrom, pos_str, snarl_id, _, ref_str_brut, alt_str_brut, p_value = fields[:7]

                # Generate placeholder sample data (e.g., "." repeated for each sample)
                sample_placeholder = "0/0"  # GT field placeholder for one sample
                placeholder = "\t".join([sample_placeholder] * num_samples)  # GT field placeholder for all samples

                # Create placeholder fields
                qual = "."
                filter_field = "PASS" if float(p_value) <= 0.05 else "LOWQ"
                info_field = f"P={p_value}"
                format_field = "GT"
                list_pos = pos_str.split(',')
                list_ref = [ref_str_brut.split(',')]
                list_alt = alt_str_brut.split(':')
                list_list_alt = [alt_str.split(',') for alt_str in list_alt]

                if snarl_id == "5262717_5262719" :
                    print("list_pos : ", list_pos)
                    print("list_ref : ", list_ref)
                    print("list_list_alt : ", list_list_alt)

                for idx, pos in enumerate(list_pos) :
                    for alt, ref in zip(list_list_alt[idx], list_ref[idx]) :
                        # Create and write the VCF line
                        vcf_line = f"{chrom}\t{pos}\t{snarl_id}\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info_field}\t{format_field}\t{placeholder}\n"
                        out_vcf.write(vcf_line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GWAS output to a VCF file using an input VCF header.")
    parser.add_argument('-g', '--gwas', type=str, required=True, help="Path to the GWAS output file (TSV format).")
    parser.add_argument('-i', '--input_vcf', type=str, required=True, help="Path to the input VCF file (to extract header information).")
    parser.add_argument('-o', '--output_vcf', type=str, required=True, help="Path to the output VCF file.")
    
    args = parser.parse_args()

    # Run the function with the parsed arguments
    create_vcf_from_gwas(args.gwas, args.input_vcf, args.output_vcf)

# python3 src/vcf_maker.py -g output/run_20241228_003954/binary_analysis.tsv -i tests/simulation/binary_data/merged_output.vcf -o test_vcf.vcf
