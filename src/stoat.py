import argparse
import list_snarl_paths
import snarl_analyser
import utils
import p_value_analysis
import write_position
import gaf_creator
import time
import logging
import os
import sys
from datetime import datetime

def main() : 

    # Argument Parsing
    parser = argparse.ArgumentParser(description='Parse and analyze snarl from VCF file')
    parser.add_argument('-p', type=utils.check_file, help='The input pangenome .pg file', required=True)
    parser.add_argument('-d', type=utils.check_file, help='The input distance index .dist file', required=True)
    parser.add_argument("-t", type=list_snarl_paths.check_threshold, help='Children threshold', required=False)
    parser.add_argument("-v", type=utils.check_format_vcf_file, help="Path to the merged VCF file (.vcf or .vcf.gz)", required=True)
    parser.add_argument("-r", type=utils.check_format_vcf_file, help="Path to the VCF file referencing all snarl positions (.vcf or .vcf.gz)", required=False)
    parser.add_argument("--listpath", type=utils.check_format_list_path, help="Path to the list paths", required=False)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--binary", type=utils.check_format_pheno, help="Path to the binary group file (.txt or .tsv)")
    group.add_argument("-q", "--quantitative", type=utils.check_format_pheno, help="Path to the quantitative phenotype file (.txt or .tsv)")
    parser.add_argument("-c", "--covariate", type=utils.check_covariate_file, required=False, help="Path to the covariate file (.txt or .tsv)")
    parser.add_argument("-g", "--gaf", required=False, help="Prepare binary gwas output to do gaf file + make gaf on 10th significatif paths")
    parser.add_argument("-o", "--output", type=str, required=False, help="Base path for the output directory")
    args = parser.parse_args()

    if args.quantitative and args.gaf:
        parser.error("The '--gaf' argument cannot be used with the '--quantitative' ('-q') argument.")

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
    list_samples = utils.parsing_samples_vcf(args.v)

    if args.covariate :
        covar = utils.parse_covariate_file(args.covariate)
        utils.check_mathing(covar, list_samples, "covariate")
    else :
        covar = None

    if args.binary:
        pheno = utils.parse_pheno_binary_file(args.binary)
        merged_dict = pheno[0].copy()  # Make a copy to avoid modifying dict1
        merged_dict.update(pheno[1]) 
        utils.check_mathing(merged_dict, list_samples, "phenotype")

    elif args.quantitative:
        pheno = utils.parse_pheno_quantitatif_file(args.quantitative)
        utils.check_mathing(pheno, list_samples, "phenotype")

    if not parser.listpath : 
        # Step 1: Parse the Pangenome Graph and Create Snarl Paths to Test
        start_time = time.time()
        logger.info("Starting snarl path decomposition...")
        stree, pg, root = list_snarl_paths.parse_graph_tree(args.p, args.d)
        snarls = list_snarl_paths.save_snarls(stree, root)
        logger.info(f"Total of snarls found : {len(snarls)}")
        logger.info("Saving snarl path decomposition...")

        output_snarl_path_not_analyse = os.path.join(output_dir, "snarl_not_analyse.tsv")
        output_snarl_path = os.path.join(output_dir, "snarl_paths.tsv")

        threshold = int(args.t) if args.t else 10 
        snarl_paths, paths_number_analysis = list_snarl_paths.loop_over_snarls_write(stree, snarls, pg, output_snarl_path, output_snarl_path_not_analyse, threshold)
        logger.info(f"Total of paths analyse : {paths_number_analysis}")

    else :
        if parser.p or parser.d : 
            logger.info("list snarls path are provided, .pg and .dist will be not analyse")
        input_snarl_path = parser.listpath
        snarl_paths, paths_number_analysis = utils.parse_snarl_path_file(input_snarl_path)
        logger.info(f"Total of snarls found : {paths_number_analysis}")

    # Step 2: Parse VCF Files and Fill the Matrix
    vcf_object = snarl_analyser.SnarlProcessor(args.v, list_samples)
    logger.info("Starting fill matrix...")
    vcf_object.fill_matrix()

    reference_vcf = args.r if args.r else args.v

    # Step 3: P-value Analysis (Binary or Quantitative)
    # Handle Binary Analysis
    if args.binary:
        gaf = True if args.gaf else False
        logger.info("Parsing binary phenotype...")
        output_snarl = os.path.join(output_dir, "binary_analysis.tsv")
        logger.info("Binary table creation...")
        vcf_object.binary_table(snarl_paths, pheno, covar, gaf, output_snarl)
        logger.info("Writing position...")
        write_position.write_pos_snarl(reference_vcf, output_snarl, "binary")

        output_manh = os.path.join(output_dir, "manhattan_plot_binary.png")
        output_qq = os.path.join(output_dir, "qq_plot_binary.png")
        output_significative = os.path.join(output_dir, "top_variant_binary.tsv")
        logger.info("Binary p-value analysis...")
        p_value_analysis.significative_snarl_binary(output_snarl, output_significative)
        p_value_analysis.qq_plot_binary(output_snarl, output_qq)
        p_value_analysis.plot_manhattan_binary(output_snarl, output_manh)
        if gaf :
            output_gaf = os.path.join(output_dir, "group_paths.gaf")
            gaf_creator.parse_input_file(output_snarl, snarl_paths, pg, output_gaf)

    # Handle Quantitative Analysis
    elif args.quantitative:
        logger.info("Parsing quantitative phenotype...")
        output_file = os.path.join(output_dir, "quantitative_analysis.tsv")
        logger.info("Quantitative table creation...")
        vcf_object.quantitative_table(snarl_paths, pheno, covar, output_file)
        logger.info("Writing position...")
        write_position.write_pos_snarl(reference_vcf, output_file, "quantitatif")

        output_manh = os.path.join(output_dir, "manhattan_plot_quantitative.png")
        output_qq = os.path.join(output_dir, "qq_plot_quantitative.png")
        output_significative = os.path.join(output_dir, "top_variant_quantitative.tsv")
        logger.info("Quantitative p-value analysis...")
        p_value_analysis.significative_snarl_quantitatif(output_file, output_significative)
        p_value_analysis.qq_plot_quantitatif(output_file, output_qq)
        p_value_analysis.plot_manhattan_quantitatif(output_file, output_manh)

    logger.info(f"GWAS analysis completed in {time.time() - start_time:.2f} seconds.")

main()

"""
Usage example:
    python3 src/stoat.py -p ../droso_data/fly/fly.pg -d ../droso_data/fly/fly.dist -v ../droso_data/pangenome.dm6.vcf \
    -r ../droso_data/fly/fly.deconstruct.vcf -q ../droso_data/pangenome_phenotype.tsv -o output

Usage test:
    python3 src/stoat.py -p test/simulation/pg.full.pg -d test/simulation/pg.dist -v test/simulation/merged_output.vcf.gz \
    -b test/simulation/phenotype.tsv -o test/simulation
"""