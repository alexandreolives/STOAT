import sys
import os
import pytest
from unittest.mock import patch
from pathlib import Path

# Add the parent directory to the Python path to enable imports from src
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

@pytest.fixture
def test_snarl_analyser():
    pg_file = "test/simulation/pg.full.pg"
    dist_file = "test/simulation/pg.dist"
    vcf_file = "test/simulation/merged_output.vcf"
    phenotype_file = "test/simulation/phenotype.tsv"
    output_dir = Path("test/simulation/simulation_output")  # Use Path for easier path handling

    # Command-line arguments to simulate
    args = [
        'src/stoat.py',
        '-p', pg_file,
        '-d', dist_file,
        '-v', vcf_file,
        '-b', phenotype_file,
        '-o', str(output_dir)
    ]

    # Mock sys.argv to simulate the command-line arguments
    with patch.object(sys, 'argv', args):
        from src.stoat import main
        main()

    # Verify that the output directory exists after execution
    assert output_dir.exists(), f"Output directory {output_dir} does not exist."

    # python3 src/stoat.py -p test/simulation/pg.full.pg -d test/simulation/pg.dist -v test/simulation/merged_output.vcf.gz -b test/simulation/phenotype.tsv -o test/simulation
