# import sys
# import os
# from pathlib import Path

# # Add the parent directory to the Python path to enable imports from src
# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

# # Import necessary modules
# import snarl_analyser
# import utils

# def test_snarl_analyser():
#     vcf_file = "tests/simulation/quantitative_data/merged_output.vcf"
#     phenotype_file = "tests/simulation/quantitative_data/phenotype.tsv"
#     snarl_file = "tests/simulation/quantitative_data/snarl_paths.tsv"
#     output_dir = Path("tests/quantitative_tests_output")

#     # Perform test logic
#     list_samples = utils.parsing_samples_vcf(vcf_file)
#     vcf_object = snarl_analyser.SnarlProcessor(vcf_file, list_samples)
#     vcf_object.fill_matrix()
#     snarl = utils.parse_snarl_path_file(snarl_file)[0]

#     os.makedirs(output_dir, exist_ok=True)
#     output = os.path.join(output_dir, "quantitative_test.assoc.tsv")

#     quantitative_pheno = utils.parse_pheno_quantitatif_file(phenotype_file)
#     vcf_object.quantitative_table(snarl, quantitative_pheno, None, False, output)

#     assert os.path.exists(output), f"Output file {output} was not created."

#     # Compare output to expected output simulation
#     expected_output = "tests/simulation/quantitative_data/expected_quantitative/quantitative_test.assoc.tsv"

#     with open(expected_output, 'r') as expected_file:
#         expected_content = expected_file.read()

#     with open(output, 'r') as output_file:
#         output_content = output_file.read()

#     assert output_content == expected_content, "The output content does not match the expected content."

# # pytest tests/test_quantitative.py
