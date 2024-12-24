import pytest
import sys
import os

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