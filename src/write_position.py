import os
import re

def get_first_snarl(s):

    match = re.findall(r'\d+', s)
    if match:
        return str(match[0])
    return None  # Return None if no integers are found

def modify_snarl(input_str):
    """
    Modifies a string by replacing '<' or '>' between numbers with an underscore.
    """
    return re.sub(r'[<>]([0-9]+)[<>]([0-9]+)', r'\1_\2', input_str)

def reverse_numbers(input_str):
    """
    Reverses the numbers on either side of an underscore in a given string.
    """
    # Split the string by the underscore
    parts = input_str.split('_')
    return f"{parts[1]}_{parts[0]}"

def classify_variant(ref, list_alt) :

    list_type_var = []
    for alt in list_alt :
        if len(ref) == len(alt) == 1:
            list_type_var.append("SNP")
        elif len(ref) > len(alt) :
            list_type_var.append("DEL")
        elif len(ref) < len(alt):
            list_type_var.append("INS")
        elif len(ref) == len(alt) and len(ref) > 1:
            list_type_var.append("MNP")
        else :
            raise ValueError(f"what is this ref : {ref}, alt : {alt}")
    
    return list_type_var
    
def write_pos_snarl(vcf_file, output_file, type):
    #vcf_dict = parse_vcf_to_dict(vcf_file)
    vcf_dict = parse_vcf_to_dict(vcf_file)
    save_info = vcf_dict.get(1, ("NA", "NA", "NA", "NA", "NA"))

    # Use a temporary file to write the updated lines
    temp_output_file = output_file + ".tmp"

    with open(output_file, 'r', encoding='utf-8') as in_f, open(temp_output_file, 'w', encoding='utf-8') as out_f:
        if type == "quantitatif" :
            out_f.write("CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tP\n")
        elif type == "binary" :
            out_f.write("CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tP_FISHER\tP_CHI2\tTOTAL_SUM\tMIN_ROW_INDEX\tNUM_COLUM\tINTER_GROUP\tAVERAGE\tGROUP_PATHS\n")
        next(in_f)
        for line in in_f:
            columns = line.strip().split('\t')
            snarl = columns[2]
            inversed_snarl = reverse_numbers(snarl)
            chrom, list_pos, list_type_var, ref, list_alt = vcf_dict.get(snarl) or vcf_dict.get(inversed_snarl) or (*save_info[:2], "NA", "NA", "UNK")

            save_info = (chrom, list_pos, list_type_var, ref, list_alt)
            columns[0], columns[4] = chrom, ref
            columns[1] = ",".join(map(str, list_pos)) if list_pos else "NA"
            columns[3] = ",".join(map(str, list_type_var)) if list_type_var != "NA" else "NA"
            columns[5] = ",".join(map(str, list_alt)) if list_alt not in ["UNK", "NA"] else list_alt

            # Write the modified line to the temp file
            out_f.write('\t'.join(columns) + '\n')

    # Replace the original file with the updated temp file
    os.replace(temp_output_file, output_file)

def write_dic(vcf_dict, fields):
    """Updates a dictionary with VCF data."""
    chr, pos, snarl_raw, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]

    # Process snarl, ref, and alt
    snarl = modify_snarl(snarl_raw)  # Pre-defined function
    ref = ref or "NA"
    alt = alt.split(",") if alt and alt != "NA" and "," in alt else [alt or "NA"]

    # Determine variant type
    variant_type = classify_variant(ref, alt) if ref != "NA" and alt != ["NA"] else ["NA"]

    # Update dictionary
    if snarl not in vcf_dict:
        vcf_dict[snarl] = [chr, [pos], variant_type, ref, alt]
    else:
        vcf_dict[snarl][1].append(pos)
        vcf_dict[snarl][2].extend(variant_type)
        vcf_dict[snarl][4] = ["NA"]  # alt becomes unknown if multiple entries exist

    return vcf_dict

def parse_vcf_to_dict(vcf_file):
    vcf_dict = {}
    number_pass = 0
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue  # Skip headers and comments

            fields = line.strip().split('\t')
            if ";" in fields[2] : # special case where ";" is find in the vcf
                decompose_snarl = fields[2].split(';')
                for snarl in decompose_snarl :
                    fields_snarl = fields[0:2] + [snarl] + fields[3:5]
                    vcf_dict = write_dic(vcf_dict, fields_snarl)
            else :
                vcf_dict = write_dic(vcf_dict, fields)

    return vcf_dict

if __name__ == "__main__" :
    reference = "/home/mbagarre/Bureau/droso_data/fly.deconstruct.vcf"
    output_file = "output/run_20241128_101551/quantitative_analysis.tsv"
    write_pos_snarl(reference, output_file)

# python3 src/write_position.py
