import argparse
import list_snarl_paths
import snarl_vcf_parser
import p_value_analysis
import write_position
import time
import logging
import os
import sys
from datetime import datetime

# Argument Parsing
parser = argparse.ArgumentParser(description='Parse and analyze snarl from VCF file')
parser.add_argument('-p', type=snarl_vcf_parser.check_file, help='The input pangenome .pg file', required=True)
parser.add_argument('-d', type=snarl_vcf_parser.check_file, help='The input distance index .dist file', required=True)
parser.add_argument("-t", type=list_snarl_paths.check_threshold, help='Children threshold', required=False)
parser.add_argument("-v", type=snarl_vcf_parser.check_format_vcf_file, help="Path to the merged VCF file (.vcf or .vcf.gz)", required=True)
parser.add_argument("-r", type=snarl_vcf_parser.check_format_vcf_file, help="Path to the VCF file referencing all snarl positions (.vcf or .vcf.gz)", required=False)

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-b", "--binary", type=snarl_vcf_parser.check_format_pheno_b, help="Path to the binary group file (.txt or .tsv)")
group.add_argument("-q", "--quantitative", type=snarl_vcf_parser.check_format_pheno_q, help="Path to the quantitative phenotype file (.txt or .tsv)")
parser.add_argument("-c", "--covariate", type=snarl_vcf_parser.check_covariate_file, required=False, help="Path to the covariate file (.txt or .tsv)")
parser.add_argument("-o", "--output", type=str, required=False, help="Base path for the output directory")
args = parser.parse_args()

# Generate unique output directory based on timestamp
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
output_dir = os.path.join(args.output or "output", f"run_{timestamp}")
os.makedirs(output_dir, exist_ok=True)

# Configure logging with both console and file handlers
log_file = os.path.join(output_dir, "run.log")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger(__name__)

# Log the exact command used to launch the script
command_line = " ".join(sys.argv)
logger.info(f"Command: {command_line}")
logger.info(f"Output directory: {output_dir}")

# Check vcf samples matching other files (pheno, covar)
list_samples = snarl_vcf_parser.parsing_samples_vcf(args.v)

if args.covariate :
    covar = snarl_vcf_parser.parse_covariate_file(args.covariate)
    snarl_vcf_parser.check_mathing(covar, list_samples, "covariate")

if args.binary:
    pheno = snarl_vcf_parser.parse_pheno_binary_file(args.binary)

elif args.quantitative:
    pheno = snarl_vcf_parser.parse_pheno_quantitatif_file(args.quantitative)

snarl_vcf_parser.check_mathing(pheno, list_samples, "phenotype")

# Step 1: Parse the Pangenome Graph and Create Snarl Paths to Test
start_time = time.time()
logger.info("Starting snarl path decomposition...")
stree, pg, root = list_snarl_paths.parse_graph_tree(args.p, args.d)
snarls = list_snarl_paths.save_snarls(stree, root)
logger.info(f"Total of snarls found : {len(snarls)}")
logger.info("Saving snarl path decomposition...")

output_snarl_not_analyse = os.path.join(output_dir, "snarl_not_analyse.tsv")

if args.t:
    snarl_paths, snarl_number_analysis = list_snarl_paths.loop_over_snarls(stree, snarls, pg, output_snarl_not_analyse, int(args.t))
else:
    snarl_paths, snarl_number_analysis = list_snarl_paths.loop_over_snarls(stree, snarls, pg, output_snarl_not_analyse)

logger.info(f"Total of paths analyse : {snarl_number_analysis}")

# Step 2: Parse VCF Files and Fill the Matrix
vcf_object = snarl_vcf_parser.SnarlProcessor(args.v, list_samples)
logger.info("Starting fill matrix...")
vcf_object.fill_matrix()

reference_vcf = args.r if args.r else args.v

# Step 3: P-value Analysis (Binary or Quantitative)
# Handle Binary Analysis
if args.binary:
    logger.info("Parsing binary phenotype...")
    output_file = os.path.join(output_dir, "binary_analysis.tsv")
    logger.info("Binary table creation...")
    vcf_object.binary_table(snarl_paths, pheno, output_file)
    logger.info("Writing position...")
    write_position.write_pos_snarl(reference_vcf, output_file)

    output_manh = os.path.join(output_dir, "manhattan_plot_binary.png")
    output_qq = os.path.join(output_dir, "qq_plot_binary.png")
    output_significative = os.path.join(output_dir, "top_variant_binary.tsv")
    logger.info("Binary p-value analysis...")
    p_value_analysis.significative_snarl_binary(output_file, output_significative)
    p_value_analysis.plot_manhattan(output_file, output_manh)
    p_value_analysis.qq_plot(output_file, output_qq)

# Handle Quantitative Analysis
elif args.quantitative:
    logger.info("Parsing quantitative phenotype...")
    output_file = os.path.join(output_dir, "quantitative_analysis.tsv")
    logger.info("Quantitative table creation...")
    vcf_object.quantitative_table(snarl_paths, pheno, output_file)
    logger.info("Writing position...")
    write_position.write_pos_snarl(reference_vcf, output_file)

    output_manh = os.path.join(output_dir, "manhattan_plot_quantitative.png")
    output_qq = os.path.join(output_dir, "qq_plot_quantitative.png")
    output_significative = os.path.join(output_dir, "top_variant_quantitative.tsv")
    logger.info("Quantitative p-value analysis...")
    p_value_analysis.significative_snarl_quantitatif(output_file, output_significative)
    p_value_analysis.plot_manhattan(output_file, output_manh)
    p_value_analysis.qq_plot(output_file, output_qq)

logger.info(f"GWAS analysis completed in {time.time() - start_time:.2f} seconds.")

"""
Usage example:
    python3 src/snarl_project.py -p ../snarl_data/fly.pg -d ../snarl_data/fly.dist -v ../droso_data/pangenome.dm6.vcf \
    -r ../droso_data/fly.normalized.vg_deconstruct.vcf -q ../droso_data/pangenome_phenotype.tsv -o output

Usage test:
    python3 src/snarl_project.py -p ...
"""