import sys
import os
import pytest
from unittest.mock import patch

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

@pytest.fixture
def mock_input_files(tmpdir):
    pg_file = tmpdir.join("test/simulation/pg.full.pg")
    dist_file = tmpdir.join("test/simulation/pg.dist")
    vcf_file = tmpdir.join("test/simulation/merged_output.vcf")
    phenotype_file = tmpdir.join("test/simulation/phenotype.tsv")
    output_dir = tmpdir.mkdir("test/simulation/simulation_output")
    
    return {
        'pg_file': pg_file,
        'dist_file': dist_file,
        'vcf_file': vcf_file,
        'phenotype_file': phenotype_file,
        'output_dir': output_dir
    }

def test_snarl_analyser(mock_input_files):
    args = [
        'src/stoat.py',
        '-p', str(mock_input_files['pg_file']),
        '-d', str(mock_input_files['dist_file']),
        '-v', str(mock_input_files['vcf_file']),
        '-b', str(mock_input_files['phenotype_file']),
        '-o', str(mock_input_files['output_dir'])
    ]

    # Mocking sys.argv to simulate the command-line arguments
    with patch.object(sys, 'argv', args):
        from src.stoat import main
        main()

    # Verify that the output directory exists
    assert mock_input_files['output_dir'].exists()

    # python3 src/stoat.py -p test/simulation/pg.full.pg -d test/simulation/pg.dist -v test/simulation/merged_output.vcf.gz -b test/simulation/phenotype.tsv -o test/simulation
