import argparse
import bdsg
from snarl_vcf_parser import parse_snarl_path_file
import re 
import math

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

def calcul_proportion_signi(number_ind_group1: int, total_ind: int, p_value: float) -> tuple:
    # Step 1: Calculate initial proportions based on a total of 60
    proportion_group1 = (number_ind_group1 / total_ind) * 60
    proportion_group2 = 60 - proportion_group1
    
    # Step 2: Calculate the adjustment factor based on a logarithmic scale of p_value
    # Adding 1e-10 to avoid log(0); adjust multiplier (e.g., 10) based on desired impact
    adjustment_factor = -math.log(max(p_value, 1e-10))  # Multiplier can be tuned for desired effect

    # Step 3: Apply the adjustment to the group with the higher initial proportion
    if proportion_group1 > proportion_group2:
        adjusted_group1 = proportion_group1 + adjustment_factor
        adjusted_group2 = proportion_group2 - adjustment_factor
    else:
        adjusted_group1 = proportion_group1 - adjustment_factor
        adjusted_group2 = proportion_group2 + adjustment_factor
    
    # Step 4: Ensure values remain within bounds [0, 60] and maintain a sum of 60
    adjusted_group1 = max(0, min(60, adjusted_group1))
    adjusted_group2 = max(0, min(60, adjusted_group2))
    
    # Rescale if needed to ensure the total sum is exactly 60
    total = adjusted_group1 + adjusted_group2
    if total != 60:
        scale_factor = 60 / total
        adjusted_group1 *= scale_factor
        adjusted_group2 *= scale_factor
    
    return int(adjusted_group1), int(adjusted_group2)

def write_gaf_lines(sequence_name : str, path : str, length : int, proportion : int, outfile : str):
    gaf_line = f"{sequence_name}\t{length}\t0\t{length}\t+\t{path}\t{length}\t0\t{length}\t{length}\t{length}\t{proportion}\tcs:Z::{length}\n"
    outfile.write(gaf_line)

def parse_input_file(input_file, snarl_dic, output_file):
    output_file_1 = add_suffix_to_filename(output_file, "_1")
    output_file_2 = add_suffix_to_filename(output_file, "_2")

    with open(input_file, 'r') as infile, open(output_file_1, 'w') as outfile1, open(output_file_2, 'w') as outfile2:
        next(infile)
        
        for line in infile:
            columns = line.strip().split()

            #CHR POS SNARL TYPE	REF	ALT	P_Fisher P_Chi2 G1_P1 G1_P2 G2_P1 G2_P2
            _, _, snarl, _, _, _, pfisher, pchi, g1p1, g1p2, g2p1, g2p2 = columns
 
            g1p1 = int(g1p1)
            g1p2 = int(g1p2)
            g2p1 = int(g2p1)
            g2p2 = int(g2p2)

            sequence_name_g1p1 = f"{snarl}_G1_{g1p1}_F{pfisher}_C{pchi}"
            sequence_name_g1p2 = f"{snarl}_G1_{g1p2}_F{pfisher}_C{pchi}"
            sequence_name_g2p1 = f"{snarl}_G2_{g2p1}_F{pfisher}_C{pchi}"
            sequence_name_g2p2 = f"{snarl}_G2_{g2p2}_F{pfisher}_C{pchi}"

            total_p1 = g1p1 + g2p1
            total_p2 = g1p2 + g2p2
            pfisher = float(pfisher)
            
            try :
                list_path = snarl_dic[snarl]
            except :
                raise ValueError(f"{snarl} not found in path list file")

            path1, path2 = list_path

            # Case when "*" is in path2
            if "*" in path2:
                less_start_path = path1
                star_path_1, star_path_2 = path2.split("*")
                star_path_1 = star_path_1[:-1]  # Remove last character from star_path_1
                length_less_start_path = calcul_path_length(pg, less_start_path)
                length_star_path_1 = calcul_path_length(pg, star_path_1)
                length_star_path_2 = calcul_path_length(pg, star_path_2)

                prop_g1_p1, prop_g2_p1 = calcul_proportion_signi(g1p1, total_p1, pfisher)
                prop_g1_p2, prop_g2_p2 = calcul_proportion_signi(g1p2, total_p2, pfisher)

                # Write lines for less_start_path
                write_gaf_lines(sequence_name_g1p1, less_start_path, length_less_start_path, prop_g1_p1, outfile1)
                write_gaf_lines(sequence_name_g2p1, less_start_path, length_less_start_path, prop_g2_p1, outfile2)

                # Write lines for star_path parts
                write_gaf_lines(sequence_name_g1p2, star_path_1, length_star_path_1, prop_g1_p2, outfile1)
                write_gaf_lines(sequence_name_g2p2, star_path_1, length_star_path_1, prop_g2_p2, outfile2)

                write_gaf_lines(sequence_name_g1p2, star_path_2, length_star_path_2, prop_g1_p2, outfile1)
                write_gaf_lines(sequence_name_g2p2, star_path_2, length_star_path_2, prop_g2_p2, outfile2)

            # Case when "*" is NOT in path2
            else:
                length_path1 = calcul_path_length(path1)
                length_path2 = calcul_path_length(path2)

                prop_g1_p1, prop_g2_p1 = calcul_proportion_signi(g1p1, total_p1, pfisher)
                prop_g1_p2, prop_g2_p2 = calcul_proportion_signi(g1p2, total_p2, pfisher)

                # Write lines for path1 for both g1p1 and g2p1
                write_gaf_lines(sequence_name_g1p1, path1, length_path1, prop_g1_p1, outfile1)
                write_gaf_lines(sequence_name_g2p1, path1, length_path1, prop_g2_p1, outfile2)

                write_gaf_lines(sequence_name_g1p2, path2, length_path2, prop_g1_p2, outfile1)
                write_gaf_lines(sequence_name_g2p2, path2, length_path2, prop_g2_p2, outfile2)
            
            break

def decompose_snarl(snarl) :
    node_list = list(map(int, re.findall(r'\d+', snarl)))
    return node_list

def calcul_path_length(pg, snarl):

    node1_id, node2_id = decompose_snarl(snarl)
    # Get the handles for the two nodes
    handle1 = pg.get_handle(node1_id)
    handle2 = pg.get_handle(node2_id)
    
    # Calculate the lengths of the nodes
    length1 = pg.get_length(handle1)
    length2 = pg.get_length(handle2)
    
    # Return the total length
    return length1 + length2

def parse_graph_tree(pg_file) :

    # load graph
    pg = bdsg.bdsg.PackedGraph()
    pg.deserialize(pg_file)

    return pg

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Parse a file and create a GAF file.")
    parser.add_argument('-s', '--snarl', type=str, help="Path to the variant snarl file.", required=True)
    parser.add_argument('-l', '--pathlist', type=str, help="Path to the list tested snarl file.", required=True)
    parser.add_argument('-p', '--pg', type=str, help='the input pangenome .pg file', required=True)
    parser.add_argument('-o', '--output', type=str, help="Path to the output GAF file.", required=False)
    args = parser.parse_args()

    pg = parse_graph_tree(args.pg)
    snarl_dic = parse_snarl_path_file(args.pathlist)

    output = args.output if args.output else 'output/output.gaf'
    parse_input_file(args.snarl, snarl_dic, output)

# python3 src/gaf_creator.py -s output/simulation_1000_binary.tsv -l ../snarl_data/simulation_1000vars_100samps/pg.snarl_netgraph.paths.tsv -p ../snarl_data/simulation_1000vars_100samps/pg.pg

