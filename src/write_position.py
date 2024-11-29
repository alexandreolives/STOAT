import os
import re
from cyvcf2 import VCF

def get_first_snarl(s):

    match = re.findall(r'\d+', s)
    if match:
        return str(match[0])
    return None  # Return None if no integers are found

def classify_variant(ref, alt) :

    if len(ref) == len(alt) == 1:
        return "SNP"
    elif len(ref) > len(alt) :
        return "DEL"
    elif len(ref) < len(alt):
        return "INS"
    elif len(ref) == len(alt) and len(ref) > 1:
        return "MNP"
    else :
        raise ValueError(f"what is this ref : {ref}, alt : {alt}")
 
def write_pos_snarl(vcf_file, output_file):
    vcf_dict = parse_vcf_to_dict(vcf_file)
    save_info = vcf_dict.get(1, ("NA", "NA", "NA", "NA", "NA"))
    
    # Use a temporary file to write the updated lines
    temp_output_file = output_file + ".tmp"

    with open(output_file, 'r', encoding='utf-8') as in_f, open(temp_output_file, 'w', encoding='utf-8') as out_f:
        out_f.write("CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tP\n")
        next(in_f)
        for line in in_f:
            columns = line.strip().split('\t')
            snarl = columns[2]
            start_snarl = snarl.split('_')[0]

            # Get VCF data or fallback to "NA" if key is not found
            chrom, pos, type_var, ref, alt = vcf_dict.get(start_snarl, (*save_info[:2], "NA", "NA", "NA"))
            save_info = (chrom, pos, type_var, ref, alt)
            columns[0], columns[1], columns[3], columns[4], columns[5] = chrom, pos, type_var, ref, alt

            # Write the modified line to the temp file
            out_f.write('\t'.join(columns) + '\n')

    # Replace the original file with the updated temp file
    os.replace(temp_output_file, output_file)

def parse_vcf_to_dict(vcf_file):
    vcf_dict = {}

    for record in VCF(vcf_file):
        # Extract VCF fields
        chr = str(record.CHROM)            # Chromosome
        pos = str(record.POS)              # Position
        snarl = get_first_snarl(record.ID) # Snarl start
        ref = str(record.REF)              # Reference allele
        alt = record.ALT[0] if record.ALT else ""
        variant_type = classify_variant(ref, alt)
        if snarl not in vcf_dict :
            vcf_dict[snarl] = (chr, pos, variant_type, ref, alt)

    return vcf_dict

if __name__ == "__main__" :
    reference = "/home/mbagarre/Bureau/droso_data/pangenome.dm6.vcf"
    output_file = "output/run_20241128_101551/quantitative_analysis.tsv"
    write_pos_snarl(reference, output_file)

# python3 src/write_position.py
