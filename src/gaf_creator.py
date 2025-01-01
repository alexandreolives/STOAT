import argparse
import bdsg # type: ignore
import utils
import re 
import math
import os 

# GAF FORMAT :
# 1     string      Query sequence name
# 2     int         Query sequence length
# 3     int         Query start (0-based; closed)
# 4     int         Query end (0-based; open)
# 5     char        Strand relative to the path: "+" or "-"
# 6     string      Path matching /([><][^\s><]+(:\d+-\d+)?)+|([^\s><]+)/
# 7     int         Path length
# 8     int         Start position on the path (0-based)
# 9     int         End position on the path (0-based)
# 10    int         Number of residue matches
# 11    int         Alignment block length
# 12    int         Mapping quality (0-255; 255 for missing)
# 13    str         unknow (cs:Z::35)

def add_suffix_to_filename(filename: str, suffix: str) -> str:
    # Split the filename into the base name and the extension
    base, ext = filename.rsplit(".", 1)
    # Add the suffix and reconstruct the filename
    new_filename = f"{base}{suffix}.{ext}"
    return new_filename

def calcul_proportion_signi(number_ind_group0: int, number_ind_group1: int, p_value: float) -> tuple:
    # Step 1: Calculate initial proportions based on a total of 60
    total_ind = number_ind_group0 + number_ind_group1
    
    if total_ind == 0 : 
        return 0,0
    
    proportion_group0 = (number_ind_group1 / total_ind) * 60
    proportion_group1 = 60 - proportion_group0
    
    # Step 2: Calculate the adjustment factor based on a logarithmic scale of p_value
    # Adding 1e-10 to avoid log(0)
    adjustment_factor = -math.log(max(p_value, 1e-10))  # Multiplier can be tuned for desired effect

    # Step 3: Apply the adjustment to the group with the higher initial proportion
    if proportion_group0 > proportion_group1:
        adjusted_group0 = proportion_group0 + adjustment_factor
        adjusted_group1 = proportion_group1 - adjustment_factor
    else:
        adjusted_group0 = proportion_group0 - adjustment_factor
        adjusted_group1 = proportion_group1 + adjustment_factor
    
    # Step 4: Ensure values remain within bounds [0, 60] and maintain a sum of 60
    adjusted_group0 = max(0, min(60, adjusted_group0))
    adjusted_group1 = max(0, min(60, adjusted_group1))
    
    # Rescale if needed to ensure the total sum is exactly 60
    total = adjusted_group0 + adjusted_group1
    if total != 60:
        scale_factor = 60 / total
        adjusted_group0 *= scale_factor
        adjusted_group1 *= scale_factor
    
    return int(adjusted_group0), int(adjusted_group1)

def write_gaf_lines(sequence_name : str, path : str, length : int, proportion : int, outfile : str):
    gaf_line = f"{sequence_name}\t{length}\t0\t{length}\t+\t{path}\t{length}\t0\t{length}\t{length}\t{length}\t{proportion}\tcs:Z::{length}\n"
    outfile.write(gaf_line)

def parse_input_file(input_file, snarl_dic, pg, output_file):
    output_file_1 = add_suffix_to_filename(output_file, "_0") # group 1
    output_file_2 = add_suffix_to_filename(output_file, "_1") # group 2

    with open(input_file, 'r') as infile, open(output_file_1, 'w') as outfile1, open(output_file_2, 'w') as outfile2:
        next(infile)
        
        for line in infile:
            columns = line.strip().split()

            #CHR POS SNARL TYPE	REF	ALT	P_Fisher P_Chi2 TOTAL_SUM MIN_ROW_INDEX NUM_COLUM INTER_GROUP AVERAGE GROUP_PATHS
            snarl = columns[2]
            pfisher, pchi = columns[6:8]
            group_paths = columns[13]
            decomposed_group_paths = group_paths.split(',')

            try :
                list_path = snarl_dic[snarl]
            except :
                raise ValueError(f"{snarl} not found in path list file")

            group_0 = []
            group_1 = []
            sequence_name_g0 = []
            sequence_name_g1 = []

            for number_group_path in decomposed_group_paths :
                path_g0, path_g1 = number_group_path.split(':')
                group_0.append(path_g0)
                group_1.append(path_g1)
                sequence_name_g0.append(f"{snarl}_G0_{path_g0}_F{pfisher}_C{pchi}")
                sequence_name_g1.append(f"{snarl}_G1_{path_g1}_F{pfisher}_C{pchi}")

            for idx, path in enumerate(list_path) :
                
                # Case where "*" is in path
                if "*" in path:
                    # Decompose the * snarl in half genereting 2 paths 
                    star_path_1, star_path_2 = path.split("*")
                    star_path_1 = star_path_1[:-1]
                    length_star_path_1 = calcul_path_length(pg, star_path_1)
                    length_star_path_2 = calcul_path_length(pg, star_path_2)
                    prop_g0, prop_g1 = calcul_proportion_signi(int(group_0[idx]), int(group_1[idx]), float(pfisher))

                    # Write lines for start_path_1
                    write_gaf_lines(sequence_name_g0[idx], star_path_1, length_star_path_1, prop_g0, outfile1)
                    write_gaf_lines(sequence_name_g1[idx], star_path_1, length_star_path_1, prop_g1, outfile2)

                    # Write lines for start_path_2
                    write_gaf_lines(sequence_name_g0[idx], star_path_2, length_star_path_2, prop_g0, outfile1)
                    write_gaf_lines(sequence_name_g1[idx], star_path_2, length_star_path_2, prop_g1, outfile2)

                # Case where "*" is NOT in path
                else:
                    length_path = calcul_path_length(pg, path)
                    prop_g0, prop_g1 = calcul_proportion_signi(int(group_0[idx]), int(group_1[idx]), float(pfisher))

                    # Write lines for path1 for both g1p1 and g2p1
                    write_gaf_lines(sequence_name_g0[idx], path, length_path, prop_g0, outfile1)
                    write_gaf_lines(sequence_name_g1[idx], path, length_path, prop_g1, outfile2)

def decompose_snarl(snarl:str) -> list:
    snarl_node = list(map(int, re.findall(r'\d+', snarl)))
    return snarl_node

def calcul_path_length(pg:bdsg.bdsg.PackedGraph, snarl:str) -> int:

    snarl_node = decompose_snarl(snarl)
    length_node = 0

    for node in snarl_node : 
        handle = pg.get_handle(node)
        length = pg.get_length(handle)
        length_node += length

    return length_node

def parse_graph_tree(pg_file:str) -> bdsg.bdsg.PackedGraph:

    # load graph
    pg = bdsg.bdsg.PackedGraph()
    pg.deserialize(pg_file)
    return pg

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Parse a file and create a GAF file.")
    parser.add_argument('-g', '--gwas', type=str, help="Path to the gwas output binary analysis file. (output file of snarl_analyser.py)", required=True)
    parser.add_argument('-l', '--pathlist', type=str, help="Path to the list tested snarl file. (output file of list_snarl_paths.py)", required=True)
    parser.add_argument('-p', '--pg', type=str, help='the input pangenome .pg file', required=True)
    parser.add_argument('-o', '--output', type=str, help="Path to the output GAF file.", required=False)
    args = parser.parse_args()

    output_dir = args.output or "output"    
    os.makedirs(output_dir, exist_ok=True)
    output = os.path.join(output_dir, "output.gaf")

    pg = parse_graph_tree(args.pg)
    snarl_dic = utils.parse_snarl_path_file(args.pathlist)[0]
    parse_input_file(args.gwas, snarl_dic, pg, output)

# python3 src/gaf_creator.py -g tests/simulation/expected_binary/binary_analysis_pos.assoc.tsv -l tests/simulation/binary_data/snarl_paths.tsv -p tests/simulation/binary_data/pg.full.pg 
