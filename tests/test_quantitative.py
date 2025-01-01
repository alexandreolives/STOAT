# import sys
# import os
# from pathlib import Path

# # Add the parent directory to the sys.path
# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

# # Import necessary modules
# import snarl_analyser
# import list_snarl_paths
# import utils

# def test_snarl_simulation_analyser():
#     vcf_file = "tests/simulation/quantitative_data/merged_output.vcf"
#     phenotype_file = "tests/simulation/quantitative_data/phenotype.tsv"
#     dist_file = "tests/simulation/quantitative_data/pg.dist"
#     pg_file = "tests/simulation/quantitative_data/pg.full.pg"
#     output_dir = Path("tests/quantitative_tests_output")

#     expected_output_snarl_path = "tests/simulation/quantitative_data/snarl_paths.tsv"
#     expected_output = "tests/simulation/quantitative_data/quantitative_analysis.assoc.tsv"

#     # Perform test logic
#     list_samples = utils.parsing_samples_vcf(vcf_file)
#     vcf_object = snarl_analyser.SnarlProcessor(vcf_file, list_samples)
#     vcf_object.fill_matrix()
#     stree, pg, root = list_snarl_paths.parse_graph_tree(pg_file,dist_file)
#     snarls = list_snarl_paths.save_snarls(stree, root)

#     output_snarl_path_not_analyse = os.path.join(output_dir, "snarl_not_analyse.tsv")
#     output_snarl_path = os.path.join(output_dir, "snarl_paths.tsv")    
#     snarl_paths = list_snarl_paths.loop_over_snarls_write(stree, snarls, pg, output_snarl_path, output_snarl_path_not_analyse, children_treshold=10)[0]

#     assert os.path.exists(output_snarl_path), f"Output file {output_snarl_path} was not created."

#     with open(expected_output_snarl_path, 'r') as expected_snarl_path:
#         expected_content_snarl_path = expected_snarl_path.read()

#     with open(output_snarl_path, 'r') as output_snarl_path:
#         output_content_output_snarl_path = output_snarl_path.read()

#     assert output_content_output_snarl_path == expected_content_snarl_path, "The output content does not match the expected content."

#     os.makedirs(output_dir, exist_ok=True)
#     output = os.path.join(output_dir, "quantitative_test.assoc.tsv")

#     quantitative_pheno = utils.parse_pheno_quantitatif_file(phenotype_file)
#     vcf_object.quantitative_table(snarl, quantitative_pheno, output=output)

#     assert os.path.exists(output), f"Output file {output} was not created."

#     with open(expected_output, 'r') as expected_file:
#         expected_content = expected_file.read()

#     with open(output, 'r') as output_file:
#         output_content = output_file.read()

#     assert output_content == expected_content, "The output content does not match the expected content."
