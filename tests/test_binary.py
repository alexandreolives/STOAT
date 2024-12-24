# import sys
# import os
# from unittest.mock import patch
# from pathlib import Path

# # Add the parent directory to the Python path to enable imports from src
# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

# def test_snarl_analyser():
#     pg_file = "tests/simulation/pg.full.pg"
#     dist_file = "tests/simulation/pg.dist"
#     vcf_file = "tests/simulation/merged_output.vcf"
#     phenotype_file = "tests/simulation/phenotype.tsv"
#     output_dir = Path("tests/simulation/simulation_output")  # Use Path for easier path handling

#     # Command-line arguments to simulate
#     args = [
#         'stoat.py',
#         '-p', pg_file,
#         '-d', dist_file,
#         '-v', vcf_file,
#         '-b', phenotype_file,
#         '-o', str(output_dir)
#     ]

#     # Mock sys.argv to simulate the command-line arguments
#     with patch.object(sys, 'argv', args):
#         from stoat import main
#         main()

#     # Verify that the output directory exists after execution
#     assert output_dir.exists(), f"Output directory {output_dir} does not exist."

# python3 tests/test_binary.py