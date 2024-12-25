import pytest
from pathlib import Path
import src.snarl_analyser
import src.utils 

@pytest.fixture
def snarl_instance():
    """Fixture to provide an instance of SnarlProcessor."""
    return src.snarl_analyser.SnarlProcessor('vcf_path', ['sample1', 'sample2'])

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
 
def test_expand_matrix(snarl_instance):
    initial_matrix = snarl_instance.matrix.get_matrix()
    initial_rows, initial_cols = initial_matrix.shape

    snarl_instance.expand_matrix()
    expanded_matrix = snarl_instance.matrix.get_matrix()
    expanded_rows, expanded_cols = expanded_matrix.shape

    assert expanded_rows == initial_rows * 2
    assert expanded_cols == initial_cols

def test_get_or_add_index(snarl_instance):
    ordered_dict = {'a': 0, 'b': 1}
    key = 'c'
    length_ordered_dict = len(ordered_dict)
    result = snarl_instance.get_or_add_index(ordered_dict, key, length_ordered_dict)
    expected = 2
    assert result == expected
    assert ordered_dict[key] == expected

    key = 'a'
    result = snarl_instance.get_or_add_index(ordered_dict, key, length_ordered_dict)
    expected = 0
    assert result == expected

def test_push_matrix(snarl_instance):
    row_header_dict = {}
    idx_snarl = 0
    decomposed_snarl = ">5123>4563"
    index_column = 0

    snarl_instance.push_matrix(idx_snarl, decomposed_snarl, row_header_dict, index_column)
    matrix = snarl_instance.matrix.get_matrix()

    assert row_header_dict[decomposed_snarl] == 0
    assert matrix[0, 0] == 1

def test_matrix(snarl_instance):

    # Use the big_vcf file from other_files directory
    vcf_path = Path('tests/other_files/big_vcf.vcf')
    list_samples = src.utils.parsing_samples_vcf(vcf_path)
    snarl_instance = src.snarl_analyser.SnarlProcessor(vcf_path, list_samples)

    # Call the method to fill the matrix
    snarl_instance.fill_matrix()

    # Verify the matrix is filled as expected
    matrix = snarl_instance.matrix.get_matrix()
    assert matrix is not None
    assert matrix.shape[0] > 0
    assert matrix.shape[1] > 0
