import sys
import os
import pytest
from unittest.mock import patch
from datetime import datetime

# Add the ../src directory to sys.path to import snarl_analyser
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import stoat
import list_snarl_paths
import snarl_analyser
import utils
import p_value_analysis
import write_position
import gaf_creator
import time
import logging

@pytest.fixture
def mock_input_files(tmpdir):
    # Create temporary files to simulate the input files
    pg_file = tmpdir.join("test/simulation/pg.full.pg")
    dist_file = tmpdir.join("test/simulation/pg.dist")
    vcf_file = tmpdir.join("test/simulation/merged_output.vcf")
    phenotype_file = tmpdir.join("test/simulation/phenotype.tsv")
    output_dir = tmpdir.mkdir("test/simulation/simulation_output")

    # Write sample content to these files
    pg_file.write("sample data")
    dist_file.write("sample data")
    vcf_file.write("sample data")
    phenotype_file.write("sample data")
    
    return {
        'pg_file': pg_file,
        'dist_file': dist_file,
        'vcf_file': vcf_file,
        'phenotype_file': phenotype_file,
        'output_dir': output_dir
    }

def test_snarl_analyser(mock_input_files):
    # Use mock input files in your test
    args = [
        'src/stoat.py',  # This is the main script you want to test
        '-p', str(mock_input_files['pg_file']),
        '-d', str(mock_input_files['dist_file']),
        '-v', str(mock_input_files['vcf_file']),
        '-b', str(mock_input_files['phenotype_file']),
        '-o', str(mock_input_files['output_dir'])
    ]

    # Mocking sys.argv to simulate the command-line arguments
    with patch.object(sys, 'argv', args):
        pass

    assert mock_input_files['output_dir'].exists()
