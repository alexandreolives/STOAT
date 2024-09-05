import unittest
import numpy as np
import snarl_vcf_parser
from unittest.mock import mock_open, patch
from collections import defaultdict

class TestSnarlProcessing(unittest.TestCase):

    def setUp(self):
        self.vcf_file = "test/small_vcf.vcf"
        self.snarl_file = 'test/list_snarl_short.txt'
        self.group_file = 'test/group.txt'
        self.pheno_file = 'test/pheno.txt'

        self.test_group_file_data = "sample\tgroup\nA\t0\nB\t1\n"
        self.test_group_file = 'test/test_group.txt'
        with open(self.test_group_file, 'w') as f:
            f.write(self.test_group_file_data)

        self.test_pheno_file_data = "FID\tIID\tPHENO\n1\tA\t0.5\n2\tB\t1.0\n"
        self.test_pheno_file = 'test/test_pheno.txt'
        with open(self.test_pheno_file, 'w') as f:
            f.write(self.test_pheno_file_data)

        self.snarl_path_file_data = "snarl\tpaths\n>1>3\t>1>2>7\n>2>3\t>1>2>3>5>6,>1>2>4>6\n"
        self.test_snarl_path_file = 'test/test_snarl_path.tsv'
        with open(self.test_snarl_path_file, 'w') as f:
            f.write(self.snarl_path_file_data)

    def test_parse_group_file(self):
        expected_group_0 = ['A']
        expected_group_1 = ['B']
        group_0, group_1 = snarl_vcf_parser.parse_group_file(self.test_group_file)
        self.assertEqual(group_0, expected_group_0)
        self.assertEqual(group_1, expected_group_1)

    def test_parse_pheno_file(self):
        expected_data = {'A': 0.5, 'B': 1.0}
        parsed_data = snarl_vcf_parser.parse_pheno_file(self.test_pheno_file)
        self.assertEqual(parsed_data, expected_data)

    def test_parse_snarl_path_file(self):
        expected_result = defaultdict(list)
        expected_result['>1>3'] = ['>1>2>7']
        expected_result['>2>3'] = ['>1>2>3>5>6', '>1>2>4>6']
        result = snarl_vcf_parser.parse_snarl_path_file(self.test_snarl_path_file)
        self.assertEqual(result, expected_result)

    def test_check_format_vcf_file(self):
        self.assertEqual(snarl_vcf_parser.check_format_vcf_file(self.vcf_file), self.vcf_file)
    
    def test_check_format_group_snarl(self):
        self.assertEqual(snarl_vcf_parser.check_format_group_snarl(self.test_group_file), self.test_group_file)
    
    def test_check_format_pheno(self):
        self.assertEqual(snarl_vcf_parser.check_format_pheno(self.test_pheno_file), self.test_pheno_file)
    
    def test_snarl_class(self):
        snarl = snarl_vcf_parser.SnarlProcessor(self.vcf_file)
        snarl.fill_matrix()
        self.assertIsInstance(snarl.matrix, snarl_vcf_parser.Matrix)

    def test_decompose_string_mixed_input(self):
        snarl = snarl_vcf_parser.SnarlProcessor(self.vcf_file)
        input_str = ">412>4<349>8"
        expected_output = [">412>4", ">4<349", "<349>8"]
        self.assertEqual(snarl.decompose_string(input_str), expected_output)

    def test_fonctionnal_binary(self):
        vcf_object = snarl_vcf_parser.SnarlProcessor(self.vcf_file)
        vcf_object.fill_matrix()
        snarls = snarl_vcf_parser.parse_snarl_path_file(self.snarl_file)
        binary_group = snarl_vcf_parser.parse_group_file(self.group_file)
        # _, list_binary_df = vcf_object.binary_table(snarl, binary_group)
        list_binary_p_value = []
        for _, list_snarl in snarls.items() :
            df = vcf_object.create_binary_table(binary_group, list_snarl)
            fisher_p_value, chi2_p_value, total_sum, inter_group, sum_column = vcf_object.binary_stat_test(df)
            list_binary_p_value.append((fisher_p_value, chi2_p_value, total_sum, inter_group, sum_column))

        expected_output = [('N/A', 'N/A', 0, 0, 0.0), (1.0, 'N/A', 1, 0, 0.5), (1.0, 'N/A', 0, 0, 0.0), ('N/A', 'N/A', 2, 1, 2.0), (1.0, 'N/A', 6, 2, 3.0), (1.0, 'N/A', 2, 1, 1.0)]
        self.assertEqual(list_binary_p_value, expected_output)

    def test_fonctionnal_pheno(self):
        vcf_object = snarl_vcf_parser.SnarlProcessor(self.vcf_file)
        vcf_object.fill_matrix()
        snarls = snarl_vcf_parser.parse_snarl_path_file(self.snarl_file)
        quantitative = snarl_vcf_parser.parse_pheno_file(self.pheno_file)
        list_quantitative_p_value= []
        for _, snarl in snarls.items() :
            df = vcf_object.create_quantitative_table(snarl)
            snarl, pvalue = vcf_object.linear_regression(df, quantitative)
            list_quantitative_p_value.append((snarl, pvalue))
        expected_output = [('>1>2>7', 'N/A'), ('>1>2>3>5>6', 0.2987424472523087), ('>1>2<2>1', 'N/A'), ('>2>3>5', 0.1173490631777916), ('>0>1>2', 0.001628646125175966), ('>1>2>3', 0.1173490631777916)]
        self.assertEqual(list_quantitative_p_value, expected_output)

if __name__ == '__main__':
    unittest.main()