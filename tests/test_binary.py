import os
from pathlib import Path

# Import necessary modules
from stoat import snarl_analyser
from stoat import list_snarl_paths
from stoat import utils
from stoat import write_position

def test_snarl_simulation_analyser():
    vcf_file = "tests/simulation/binary_data/merged_output.vcf"
    phenotype_file = "tests/simulation/binary_data/phenotype.tsv"
    dist_file = "tests/simulation/binary_data/pg.dist"
    pg_file = "tests/simulation/binary_data/pg.full.pg"
    output_dir = Path("tests/binary_tests_output")

    expected_output_snarl_path = "tests/simulation/binary_data/snarl_paths.tsv"
    expected_output = "tests/simulation/binary_data/binary_analysis.assoc.tsv"

    # Perform test logic
    list_samples = utils.parsing_samples_vcf(vcf_file)
    vcf_object = snarl_analyser.SnarlProcessor(vcf_file, list_samples)
    vcf_object.fill_matrix()
    stree, pg, root = list_snarl_paths.parse_graph_tree(pg_file,dist_file)
    snarls = list_snarl_paths.save_snarls(stree, root)

    output_snarl_path_not_analyse = os.path.join(output_dir, "snarl_not_analyse.tsv")
    output_snarl_path = os.path.join(output_dir, "snarl_paths.tsv")    
    snarl_paths = list_snarl_paths.loop_over_snarls_write(stree, snarls, pg, output_snarl_path, output_snarl_path_not_analyse, children_treshold=10)[0]

    assert os.path.exists(output_snarl_path), f"Output file {output_snarl_path} was not created."

    with open(expected_output_snarl_path, 'r') as expected_snarl_path:
        expected_content_snarl_path = expected_snarl_path.read()

    with open(output_snarl_path, 'r') as output_snarl_path:
        output_content_output_snarl_path = output_snarl_path.read()

    assert output_content_output_snarl_path == expected_content_snarl_path, "The output content does not match the expected content."

    os.makedirs(output_dir, exist_ok=True)
    output = os.path.join(output_dir, "binary_test.assoc.tsv")

    binary_group = utils.parse_pheno_binary_file(phenotype_file)
    vcf_object.binary_table(snarl_paths, binary_group, gaf=True, output=output)
    write_position.write_pos_snarl(vcf_file, output, "binary")

    assert os.path.exists(output), f"Output file {output} was not created."

    with open(expected_output, 'r') as expected_file:
        expected_content = expected_file.read()

    with open(output, 'r') as output_file:
        output_content = output_file.read()

    assert output_content == expected_content, "The output content does not match the expected content."

# pytest tests/test_binary.py
