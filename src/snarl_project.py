import argparse
import list_snarl_paths
import snarl_vcf_parser
import p_value_analysis
import time
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Argument Parsing
parser = argparse.ArgumentParser('List path through the netgraph of each snarl in a pangenome')
parser.add_argument('-p', help='The input pangenome .pg file', required=True)
parser.add_argument('-d', help='The input distance index .dist file', required=True)
parser.add_argument("-t", type=list_snarl_paths.check_threshold, help='Children threshold', required=False)
parser.add_argument("-v", type=snarl_vcf_parser.check_format_vcf_file, help="Path to the merged VCF file (.vcf or .vcf.gz)", required=True)
parser.add_argument("-r", type=snarl_vcf_parser.check_format_vcf_file, help="Path to the VCF file referencing all snarl positions (.vcf or .vcf.gz)")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-b", "--binary", type=snarl_vcf_parser.check_format_pheno_b, help="Path to the binary group file (.txt or .tsv)")
group.add_argument("-q", "--quantitative", type=snarl_vcf_parser.check_format_pheno_q, help="Path to the quantitative phenotype file (.txt or .tsv)")
parser.add_argument("-o", "--output", type=str, required=False, help="Path to the output file")
args = parser.parse_args()

# Main Process
logger.info("Starting snarl path analysis...")

# Step 1: Parse the Pangenome Graph and Create Snarl Paths
start_time = time.time()
stree, pg, root = list_snarl_paths.parse_graph_tree(args.p, args.d)
snarls = list_snarl_paths.create_snarls(stree, root)
logger.info(f"Parsed graph and created snarls in {time.time() - start_time:.2f} seconds.")

if args.t:
    snarl_paths = list_snarl_paths.loop_over_snarls(stree, snarls, pg, args.t)
else:
    snarl_paths = list_snarl_paths.loop_over_snarls(stree, snarls, pg)

logger.info(f"Generated snarl paths in {time.time() - start_time:.2f} seconds.")

# Step 2: Parse VCF Files and Fill the Matrix
start_time = time.time()
vcf_dict = snarl_vcf_parser.parse_vcf_to_dict(args.r)
vcf_object = snarl_vcf_parser.SnarlProcessor(args.v, vcf_dict)
vcf_object.fill_matrix()
logger.info(f"Matrix populated in {time.time() - start_time:.2f} seconds.")

# Step 3: P-value Analysis (Binary or Quantitative)
start_time = time.time()

if args.binary:
    binary_group = snarl_vcf_parser.parse_group_file(args.binary)
    output_file = args.output or "output/binary_analysis.tsv"
    vcf_object.binary_table(snarl_paths, binary_group, output_file)
    snarl_vcf_parser.write_pos_snarl(args.vcf_pangenome, output_file)
    logger.info("Binary analysis table created.")

    output_manh = "output/manhattan_plot_binary.png"
    output_qq = "output/qq_plot_binary.png"
    output_significative = "output/top_variant_binary.tsv"
    p_value_analysis.significative_snarl_binary(output_file, output_significative)
    p_value_analysis.plot_manhattan(output_file, output_manh)
    p_value_analysis.qq_plot(output_file, output_qq)
    logger.info("Binary p-value analysis completed with Manhattan and QQ plots.")

if args.quantitative:
    quantitative_pheno = snarl_vcf_parser.parse_pheno_file(args.quantitative)
    output_file = args.output or "output/quantitative_analysis.tsv"
    vcf_object.quantitative_table(snarl_paths, quantitative_pheno, output_file)
    snarl_vcf_parser.write_pos_snarl(args.r, output_file)
    logger.info("Quantitative analysis table created.")

    output_manh = "output/manhattan_plot_quantitative.png"
    output_qq = "output/qq_plot_quantitative.png"
    output_significative = "output/top_variant_quantitative.tsv"
    p_value_analysis.significative_snarl_quantitatif(output_file, output_significative)
    p_value_analysis.plot_manhattan(output_file, output_manh)
    p_value_analysis.qq_plot(output_file, output_qq)
    logger.info("Quantitative p-value analysis completed with Manhattan and QQ plots.")

logger.info(f"P-value analysis completed in {time.time() - start_time:.2f} seconds.")

"""
    python3 src/snarl_project.py -p ../../snarl_data/fly.pg -d ../../snarl_data/fly.dist -v ../../snarl_data/fly.merged.vcf \
    -r ../../snarl_data/fly.vcf -q ../../snarl_data/updated_phenotypes.txt -o output/simulation_quantitative.tsv
"""
