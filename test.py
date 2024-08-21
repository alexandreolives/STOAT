import unittest
import numpy as np
import pandas as pd
from unittest.mock import patch
from snarl_matrix import Chunk, Matrix, Snarl

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
        chunk3 = Chunk(3, 2)
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
    @patch('test/test_variant.vcf')
    def test_is_valid_vcf(self):
        snarl = Snarl("test/test_variant.vcf", "test/test_group.txt", "test/test_path.txt")
        self.assertTrue(snarl.is_valid_vcf("test/test_variant.vcf"))
        self.assertFalse(snarl.is_valid_vcf("test.txt"))

    def test_determine_str(self):
        snarl = Snarl("test/test_variant.vcf", "test/test_group.txt", "test/test_path.txt")
        i, substring = snarl.determine_str("1234>5678", 9, 4)
        self.assertEqual(i, 9)
        self.assertEqual(substring, "5678")

    def test_decompose_string(self):
        snarl = Snarl("test/test_variant.vcf", "test/test_group.txt", "test/test_path.txt")
        result = snarl.decompose_string(">1272>1273>1274")
        expected = [">1272>1273", ">1273>1274"]
        self.assertEqual(result, expected)

    def test_decompose_snarl(self):
        snarl = Snarl("test/test_variant.vcf", "test/test_group.txt", "test/test_path.txt")
        result = snarl.decompose_snarl([">1272>1273>1274", ">1272>1274"])
        expected = [[">1272>1273", ">1273>1274"], [">1272>1274"]]
        self.assertEqual(result, expected)

    @patch('snarl_matrix.Snarl._parse_group_file')
    @patch('snarl_matrix.Snarl._parse_path_file')
    @patch('snarl_matrix.Snarl._parse_vcf')
    def test_snarl_initialization(self, mock_parse_vcf, mock_parse_path_file, mock_parse_group_file):
        mock_parse_vcf.return_value = (np.zeros((2, 2)), [], [], [])
        mock_parse_path_file.return_value = []
        mock_parse_group_file.return_value = ([], [])

        snarl = Snarl("test/test_variant.vcf", "test/test_group.txt", "test/test_path.txt")
        self.assertIsInstance(snarl.matrix, Matrix)

    @patch('snarl_matrix.Snarl.create_table')
    def test_create_tables(self, mock_create_table):
        snarl = Snarl("test/test_variant.vcf", "test/test_group.txt", "test/test_path.txt")
        mock_create_table.return_value = pd.DataFrame({"G0": [1, 2], "G1": [3, 4]})
        result = snarl.create_tables()
        self.assertEqual(len(result), len(snarl.snarl_paths))
        self.assertIsInstance(result[0], pd.DataFrame)

class TestFonctionnal(unittest.TestCase):
    def test_create_tables(self, mock_create_table):
        snarl = Snarl("test/test_variant.vcf", "test/test_group.txt", "test/test_path.txt")

if __name__ == '__main__':
    unittest.main()