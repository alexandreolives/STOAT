import sys
import os
from unittest.mock import patch
from pathlib import Path

# Add the parent directory to the Python path to enable imports from src
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

def test_snarl_analyser():
    vcf_file = "tests/simulation/merged_output.vcf"
    phenotype_file = "tests/simulation/phenotype.tsv"
    snarl_file = "tests/simulation/snarl_paths.tsv"
    output_dir = Path("tests/simulation/simulation_output")

    # Command-line arguments to simulate
    args = [
        'snarl_analyser.py',
        vcf_file,
        snarl_file,
        '-b', phenotype_file,
        '-o', str(output_dir)
    ]

    # Mock sys.argv to simulate the command-line arguments
    with patch.object(sys, 'argv', args):
        import snarl_analyser
        snarl_analyser.main()

    assert output_dir.exists(), f"Output directory {output_dir} does not exist."

# TODO : Compare output to expected output simulation 

# pytest tests/test_binary.py
