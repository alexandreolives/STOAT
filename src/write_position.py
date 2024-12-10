import os
import re
from cyvcf2 import VCF

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
    #vcf_dict = parse_vcf_to_dict(vcf_file)
    vcf_dict = parse_vcf_to_dict_special(vcf_file)
    save_info = vcf_dict.get(1, ("NA", "NA", "NA", "NA", "NA"))

    # Use a temporary file to write the updated lines
    temp_output_file = output_file + ".tmp"

    with open(output_file, 'r', encoding='utf-8') as in_f, open(temp_output_file, 'w', encoding='utf-8') as out_f:
        out_f.write("CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tP\n")
        next(in_f)
        for line in in_f:
            columns = line.strip().split('\t')
            snarl = columns[2]
            inversed_snarl = reverse_numbers(snarl)
            chrom, list_pos, list_type_var, list_ref, list_alt = vcf_dict.get(snarl) or vcf_dict.get(inversed_snarl) or (*save_info[:2], "NA", "NA", "NA")

            # Get VCF data or fallback to "NA" if key is not found
            #chrom, list_pos, type_var, ref, alt = vcf_dict.get(snarl, (*save_info[:2], "NA", "NA", "NA"))
            save_info = (chrom, list_pos, list_type_var, list_ref, list_alt)
            #columns[0], columns[3], columns[4], columns[5] = chrom, type_var, ref, alt
            columns[0] = chrom
            columns[1] = ""
            columns[3] = ""
            columns[4] = ""
            columns[5] = ""

            # if len(list_pos) > 1 :
            #     print("save_info : ", save_info)

            for i in range(len(list_pos)):
                columns[1] += str(list_pos[i]) + ","  # Concatenate chromosome strings with commas
                try :
                    columns[3] += str(list_type_var[i]) + ","  # Concatenate chromosome strings with commas
                except IndexError :
                    print(snarl)
                    print(list_pos)
                    print(list_type_var)
                    print(list_ref)
                    exit()

                columns[4] += str(list_ref[i]) + ","  # Concatenate chromosome strings with commas
                columns[5] += str(list_alt[i]) + ","  # Concatenate chromosome strings with commas

            columns[1] = columns[1][:-1]  # Remove the trailing comma
            columns[3] = columns[3][:-1]  # Remove the trailing comma
            columns[4] = columns[4][:-1]  # Remove the trailing comma
            columns[5] = columns[5][:-1]  # Remove the trailing comma

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

def write_dic(vcf_dict, fields) :

    chr = fields[0]          # Chromosome
    pos = fields[1]          # Position
    snarl = modify_snarl(fields[2])  # Use pre-defined function
    ref = fields[3] if fields[3] is not None and fields[3] != "" else "NA"
    alt = fields[4] if fields[4] is not None and fields[4] != "" else "NA"
    if alt != 'NA' and len(alt) > 1 :
        alt = alt.split(",")[0]
        
    if ref != 'NA' and alt != 'NA' :
        variant_type = classify_variant(ref, alt)  # Use pre-defined function
    else :
        variant_type = "NA"

    if snarl not in vcf_dict:
        vcf_dict[snarl] = [chr, [pos], [variant_type], [ref], [alt]]

    else :
        vcf_dict[snarl][1].append(pos)
        vcf_dict[snarl][2].append(variant_type)
        vcf_dict[snarl][3].append(ref)
        vcf_dict[snarl][4].append(alt)
    return vcf_dict

def parse_vcf_to_dict_special(vcf_file):
    vcf_dict = {}
    number_pass = 0
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue  # Skip headers and comments

            fields = line.strip().split('\t')
            if ";" in fields[2] : # special case where ";" is find in the vcf
                decompose_snarl = fields[2].split(';')
                fields_1 = fields[0:2] + [decompose_snarl[0]] + fields[3:5]
                fields_2 = fields[0:2] + [decompose_snarl[1]] + fields[3:5]
                vcf_dict = write_dic(vcf_dict, fields_1)
                vcf_dict = write_dic(vcf_dict, fields_2)
            else :
                vcf_dict = write_dic(vcf_dict, fields)

    print("number of variant not parsed : ", number_pass)
    return vcf_dict

if __name__ == "__main__" :
    reference = "/home/mbagarre/Bureau/pangenome.dm6.normalized.vcf"
    output_file = "output/run_20241128_101551/quantitative_analysis.tsv"
    write_pos_snarl(reference, output_file)

# python3 src/write_position.py
