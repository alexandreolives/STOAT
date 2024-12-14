import unittest
from typing import List
import sys
import os

# Add the ../src directory to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import snarl_analyser

class TestDetermineStr(unittest.TestCase):
    
    def test_decompose_string(self):
        # Instance of the class containing the function
        SnarlClass = snarl_analyser.SnarlProcessor('tt','t')
        
        # Test case 1
        s = ">5123>4563<23789"
        result = SnarlClass.decompose_string(s)
        expected = [">5123>4563", ">4563<23789"]
        self.assertEqual(result, expected)

        # Test case 2
        s = "<1>522<3335"
        result = SnarlClass.decompose_string(s)
        expected = ["<1>522", ">522<3335"]
        self.assertEqual(result, expected)

    def test_decompose_snarl(self):
        # Instance of the class containing the function
        SnarlClass = snarl_analyser.SnarlProcessor('tt','t')
        
        # Test case 1
        lst = [">5123>4563<23789", ">1>522<333"]
        result = SnarlClass.decompose_snarl(lst)
        expected = [
            [">5123>4563", ">4563<23789"],
            [">1>522", ">522<333"]
        ]
        self.assertEqual(result, expected)

        # Test case 2: List with empty strings
        lst = ["", "<5123>4563"]
        result = SnarlClass.decompose_snarl(lst)
        expected = [
            [],
            ["<5123>4563"]
        ]
        self.assertEqual(result, expected)

if __name__ == "__main__":
    unittest.main()
