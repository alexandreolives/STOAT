import unittest
import numpy as np
import pandas as pd
from unittest.mock import patch
import snarl_vcf_parser
from snarl_vcf_parser import Chunk, Matrix, Snarl
import argparse

class TestChunk(unittest.TestCase):
    def test_add_data(self):
        chunk = Chunk(2, 2)
        chunk.add_data(0, 0)
        self.assertTrue(chunk.get_data()[0, 0])
        self.assertFalse(chunk.get_data()[1, 1])

    def test_concatenate_matrix(self):
        chunk1 = Chunk(2, 2)
        chunk1.add_data(0, 0)
        chunk2 = Chunk(2, 2)
        chunk2.add_data(1, 1)

        concatenated = Chunk.concatenate_matrix([chunk1, chunk2], axis=0)
        self.assertEqual(concatenated.shape, (4, 2))
        self.assertTrue(concatenated[0, 0])
        self.assertTrue(concatenated[3, 1])

        # Test concatenation when input list is empty
        with self.assertRaises(ValueError):
            Chunk.concatenate_matrix([], axis=0)
        
        # Test when shapes are incompatible
        chunk3 = Chunk(3, 3)
        with self.assertRaises(ValueError):
            Chunk.concatenate_matrix([chunk1, chunk3], axis=0)

class TestMatrix(unittest.TestCase):
    def test_matrix_initialization(self):
        matrix = np.zeros((2, 2))
        row_header = ["row1", "row2"]
        column_header = ["col1", "col2"]
        list_samples = ["sample1", "sample2"]

        matrix_obj = Matrix((matrix, row_header, column_header, list_samples))
        self.assertTrue(np.array_equal(matrix_obj.get_matrix(), matrix))
        self.assertEqual(matrix_obj.get_row_header(), row_header)
        self.assertEqual(matrix_obj.get_column_header(), column_header)
        self.assertEqual(matrix_obj.get_list_samples(), list_samples)

class TestSnarl(unittest.TestCase):
    def test_is_valid_vcf(self):
        self.assertTrue(snarl_vcf_parser.check_format_vcf_file("../test/test_variant.vcf"))
        
        with self.assertRaises(argparse.ArgumentTypeError):
            self.assertFalse(snarl_vcf_parser.check_format_vcf_file("../test/test.txt"))

    def test_decompose_snarl(self):
        snarl = Snarl("../test/test_variant.vcf")
        result = snarl.decompose_snarl([">1322<1323>1323<1323>1325", ">1272>1274"])
        expected = [['>1322<1323', '<1323>1323', '>1323<1323', '<1323>1325'], [">1272>1274"]]
        self.assertEqual(result, expected)

    def test_decompose_star_snarl(self):
        snarl = Snarl("../test/test_variant.vcf")
        result = snarl.decompose_snarl([">1272<1273>*>1276<1277", ">1272>1274"])
        expected = [[">1272<1273", "<1273>*", ">*>1276", ">1276<1277"], [">1272>1274"]]
        self.assertEqual(result, expected)

    @patch('snarl_vcf_parser.Snarl.create_matrix')
    def test_snarl_initialization(self, mockcreate_matrix):
        mockcreate_matrix.return_value = (np.zeros((2, 2)), [], [], [])
        snarl = Snarl("../test/test_variant.vcf")
        self.assertEqual(snarl.matrix, None)

class TestFonctionnal(unittest.TestCase) :
    def test_binary(self) :
        vcf_object = Snarl("../test/test_variant.vcf")
        vcf_object.initialise_matrix()
        snarl = snarl_vcf_parser.parse_snarl_path_file("../test/test_path.txt")
        binary_group = snarl_vcf_parser.parse_group_file("../test/test_group.txt")
        binary_tables = vcf_object.binary_table(snarl, binary_group)
        binary_tables_expected = pd.DataFrame()

if __name__ == '__main__':
    unittest.main()