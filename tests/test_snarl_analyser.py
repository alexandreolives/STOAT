import pytest
import sys
import os
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import snarl_analyser

@pytest.fixture
def snarl_instance():
    """Fixture to provide an instance of SnarlProcessor."""
    return snarl_analyser.SnarlProcessor('vcf_path', ['sample1', 'sample2'])

def test_decompose_string(snarl_instance):
    # Test case 1
    s = ">5123>4563<23789"
    result = snarl_instance.decompose_string(s)
    expected = [">5123>4563", ">4563<23789"]
    assert result == expected

    # Test case 2
    s = "<1>522<3335"
    result = snarl_instance.decompose_string(s)
    expected = ["<1>522", ">522<3335"]
    assert result == expected

def test_decompose_snarl(snarl_instance):
    # Test case 1
    lst = [">5123>4563<23789", ">1>522<333"]
    result = snarl_instance.decompose_snarl(lst)
    expected = [
        [">5123>4563", ">4563<23789"],
        [">1>522", ">522<333"]
    ]
    assert result == expected

    # Test case 2: List with empty strings
    lst = ["", "<5123>4563"]
    result = snarl_instance.decompose_snarl(lst)
    expected = [
        [],
        ["<5123>4563"]
    ]
    assert result == expected

def test_quantitative_table(snarl_instance):
    # Mock data for snarls and quantitative
    snarls = {
        'snarl1': ['>5123>4563', '>4563<23789'],
        'snarl2': ['<1>522', '>522<3335']
    }
    quantitative = {
        'sample1': 1.5,
        'sample2': 2.5
    }

    # Call the method
    snarl_instance.matrix.set_matrix(pd.DataFrame([[1, 0], [0, 1]], columns=['snarl1', 'snarl2']))
    snarl_instance.matrix.set_row_header({'snarl1': 0, 'snarl2': 1})
    snarl_instance.list_samples = ['sample1', 'sample2']
    snarl_instance.vcf_path = 'vcf_path'
    snarl_instance.quantitative_table(snarls, quantitative, output="output/quantitative_output.tsv")

    # Read the output file and check the content
    with open("output/quantitative_output.tsv", 'r') as f:
        lines = f.readlines()
        assert lines[0].strip() == 'CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tRSQUARED\tBETA\tSE\tP'
        assert len(lines) > 1  # Ensure there are results written

    snarl_instance.binary_table(snarls, quantitative, output="output/quantitative_output.tsv")

def test_binary_table(snarl_instance):
    # Mock data for snarls and binary groups
    snarls = {
        'snarl1': ['>5123>4563', '>4563<23789'],
        'snarl2': ['<1>522', '>522<3335']
    }
    binary_groups = {
        'sample1' : 0,
        'sample2' : 1
    }

    # Call the method
    snarl_instance.matrix.set_matrix(pd.DataFrame([[1, 0], [0, 1]], columns=['snarl1', 'snarl2']))
    snarl_instance.matrix.set_row_header({'snarl1': 0, 'snarl2': 1})
    snarl_instance.list_samples = ['sample1', 'sample2']
    snarl_instance.vcf_path = 'vcf_path'
    snarl_instance.binary_table(snarls, binary_groups, output="output/binary_output.tsv")

    # Read the output file and check the content
    with open("output/binary_output.tsv", 'r') as f:
        lines = f.readlines()
        assert lines[0].strip() == 'CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tP_FISHER\tP_CHI2\tTOTAL_SUM\tMIN_ROW_INDEX\tNUM_COLUM\tINTER_GROUP\tAVERAGE\tGROUP_PATHS'
        assert len(lines) > 1  # Ensure there are results written

# python3 tests/test_snarl_analyser.py
