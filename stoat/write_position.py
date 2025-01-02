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

    vcf_dict = parse_vcf_to_dict(vcf_file)
    save_info = vcf_dict.get(1, ("NA", {"NA" : {"NA" : "NA"}}, "NA"))
    seen = False

    # Use a temporary file to write the updated lines
    temp_output_file = output_file + ".tmp"

    with open(output_file, 'r', encoding='utf-8') as in_f, open(temp_output_file, 'w', encoding='utf-8') as out_f:
        if type == "quantitatif" :
            out_f.write("CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tRSQUARED\tBETA\tSE\tP\n")
        elif type == "binary" :
            out_f.write("CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tP_FISHER\tP_CHI2\tTOTAL_SUM\tMIN_ROW_INDEX\tNUM_COLUM\tINTER_GROUP\tAVERAGE\tGROUP_PATHS\n")
        next(in_f)
        for line in in_f:
            columns = line.strip().split('\t')
            snarl = columns[2]
            inversed_snarl = reverse_numbers(snarl)
            chrom, dict_pos_ref_alt, list_type_var = vcf_dict.get(snarl) or vcf_dict.get(inversed_snarl) or (save_info[0], save_info[1], save_info[2])

            pos_dic = {}
            for pos, ref_alt_dict in dict_pos_ref_alt.items():
                for ref_base, alt_list in ref_alt_dict.items():
                    if not pos in pos_dic :
                        pos_dic[pos] = {ref_base : alt_list}
                    else :
                        pos_dic[pos].update({ref_base : alt_list})

            save_info = (chrom, dict_pos_ref_alt, list_type_var)
            columns[0] = chrom
            columns[1] = ",".join(map(str, pos_dic.keys()))  # POS
            columns[3] = ",".join(map(str, list_type_var))   # TYPE VAR (INS, DEL, etc.)
            columns[4] = ":".join(",".join(map(str, pos_dic[pos].keys())) for pos in pos_dic)  # REF
            # ALT: list of lists formatted with ',' for mutltiple alt in 1 ref and ':' for each position
            columns[5] = ":".join([",".join(map(str, alt_list)) for pos in pos_dic for alt_list in pos_dic[pos].values()])

            # Write the modified line to the temp file
            out_f.write('\t'.join(columns) + '\n')

    # Replace the original file with the updated temp file
    os.replace(temp_output_file, output_file)

def write_dic(vcf_dict, fields):
    """Updates a dictionary with VCF data."""
    chr, pos, snarl_raw, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]

    # Process snarl, ref, and alt
    snarl = modify_snarl(snarl_raw)  # Pre-defined function
    alt = alt.split(",") if "," in alt else [alt]

    # Determine variant type
    variant_type = classify_variant(ref, alt)

    # Update dictionary
    if snarl not in vcf_dict:
        dict_pos_ref_alt = {pos : {ref : alt}}
        vcf_dict[snarl] = [chr, dict_pos_ref_alt, variant_type]
    else:
        if pos not in vcf_dict[snarl][1]:
            vcf_dict[snarl][1][pos] = {ref: alt}
        else:
            if ref not in vcf_dict[snarl][1][pos]:
                vcf_dict[snarl][1][pos][ref] = alt
            else:
                vcf_dict[snarl][1][pos][ref].extend(alt)
        vcf_dict[snarl][2].extend(variant_type)

    return vcf_dict

def parse_vcf_to_dict(vcf_file):
    vcf_dict = {}
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
    reference = "../droso_data/pangenome.dm6.vcf"
    output_file = "output/run_20250102_134046/quantitative_analysis.tsv"
    write_pos_snarl(reference, output_file, "quantitatif")

# python3 stoat/write_position.py
