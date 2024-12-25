import os
from pathlib import Path

# Import necessary modules
import src.snarl_analyser
import src.utils

def test_snarl_analyser():
    vcf_file = "tests/simulation/binary_data/merged_output.vcf"
    phenotype_file = "tests/simulation/binary_data/phenotype.tsv"
    snarl_file = "tests/simulation/binary_data/snarl_paths.tsv"
    output_dir = Path("tests/binary_tests_output")
    expected_output = "tests/simulation/expected_binary/binary_analysis.assoc.tsv"

    # Perform test logic
    list_samples = src.utils.parsing_samples_vcf(vcf_file)
    vcf_object = src.snarl_analyser.SnarlProcessor(vcf_file, list_samples)
    vcf_object.fill_matrix()
    snarl = src.utils.parse_snarl_path_file(snarl_file)[0]

    os.makedirs(output_dir, exist_ok=True)
    output = os.path.join(output_dir, "binary_test.assoc.tsv")

    binary_group = src.utils.parse_pheno_binary_file(phenotype_file)
    vcf_object.binary_table(snarl, binary_group, None, True, output)

    assert os.path.exists(output), f"Output file {output} was not created."

    with open(expected_output, 'r') as expected_file:
        expected_content = expected_file.read()

    with open(output, 'r') as output_file:
        output_content = output_file.read()

    assert output_content == expected_content, "The output content does not match the expected content."

# pytest tests/test_binary.py
