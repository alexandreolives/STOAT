import unittest
import numpy as np
import pandas as pd
from unittest.mock import patch
import snarl_vcf_parser
from collections import defaultdict

class TestSnarlProcessing(unittest.TestCase):

    def setUp(self):
        self.vcf_file = "test/small_vcf.vcf"

        self.group_file_data = "A\tB\nC\tD\n"
        self.group_file = 'test_group.txt'
        with open(self.group_file, 'w') as f:
            f.write(self.group_file_data)

        self.pheno_file_data = "FID\tIID\tPHENO\n1\tA\t0.5\n2\tB\t1.0\n"
        self.pheno_file = 'test_pheno.txt'
        with open(self.pheno_file, 'w') as f:
            f.write(self.pheno_file_data)

        self.snarl_path_file_data = "snarl\tpaths\n>1>3\t>1>2>7\n>2>3\t>1>2>3>5>6,>1>2>4>6\n"
        self.snarl_path_file = 'test_snarl_path.tsv'
        with open(self.snarl_path_file, 'w') as f:
            f.write(self.snarl_path_file_data)

    def test_parse_group_file(self):
        expected_group_0 = ['A', 'C']
        expected_group_1 = ['B', 'D']
        group_0, group_1 = snarl_vcf_parser.parse_group_file(self.group_file)
        self.assertEqual(group_0, expected_group_0)
        self.assertEqual(group_1, expected_group_1)

    def test_parse_pheno_file(self):
        expected_data = {'A': 0.5, 'B': 1.0}
        parsed_data = snarl_vcf_parser.parse_pheno_file(self.pheno_file)
        self.assertEqual(parsed_data, expected_data)

    def test_parse_snarl_path_file(self):
        expected_result = defaultdict(list)
        expected_result['>1>3'] = ['>1>2>7']
        expected_result['>2>3'] = ['>1>2>3>5>6', '>1>2>4>6']
        result = snarl_vcf_parser.parse_snarl_path_file(self.snarl_path_file)
        self.assertEqual(result, expected_result)

    def test_check_format_vcf_file(self):
        self.assertEqual(snarl_vcf_parser.check_format_vcf_file(self.vcf_file), self.vcf_file)
    
    def test_check_format_group_snarl(self):
        self.assertEqual(snarl_vcf_parser.check_format_group_snarl(self.group_file), self.group_file)
    
    def test_check_format_pheno(self):
        self.assertEqual(snarl_vcf_parser.check_format_pheno(self.pheno_file), self.pheno_file)
    
    def test_chunk(self):
        chunk = snarl_vcf_parser.Chunk(MATRIX_ROW_NUMBER=2, nb_matrix_column=3)
        chunk.add_data(0, 1)
        expected_chunk = np.array([[False, True, False], [False, False, False]])
        np.testing.assert_array_equal(chunk.get_chunk(), expected_chunk)

    def test_matrix(self):
        matrix_data = np.array([[True, False], [False, True]])
        row_header = ['row1', 'row2']
        column_header = ['col1', 'col2']
        matrix = snarl_vcf_parser.Matrix(matrix_data, row_header, column_header)
        self.assertTrue(np.array_equal(matrix.get_matrix(), matrix_data))
        self.assertEqual(matrix.get_row_header(), row_header)
        self.assertEqual(matrix.get_column_header(), column_header)

    def test_snarl_class(self):
        snarl = snarl_vcf_parser.SnarlProcessor(self.vcf_file)
        snarl.fill_matrix()
        self.assertIsInstance(snarl.matrix, snarl_vcf_parser.Matrix)

    def test_decompose_string_mixed_input(self):
        snarl = snarl_vcf_parser.SnarlProcessor(self.vcf_file)
        input_str = ">412>4<349>8"
        expected_output = [">412>4", ">4<349", "<349>8"]
        self.assertEqual(snarl.decompose_string(input_str), expected_output)

if __name__ == '__main__':
    unittest.main()