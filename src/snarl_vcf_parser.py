import argparse
from cyvcf2 import VCF
from typing import List
import numpy as np
import pandas as pd
import statsmodels.api as sm
from collections import defaultdict
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
import os
import re
import time

class Matrix :
    def __init__(self, default_row_number=1_000_000, column_number=2):
        self.default_row_number = default_row_number 
        self.matrix = np.zeros((default_row_number, column_number),dtype=bool)
        self.row_header = None
    
    def get_matrix(self):
        return self.matrix

    def set_matrix(self, expended_matrix) :
        self.matrix = expended_matrix

    def get_row_header(self):
        return self.row_header

    def get_default_row_number(self):
        return self.default_row_number

    def set_row_header(self, row_header):
        self.row_header = row_header  

    def add_data(self, idx_snarl, idx_geno):
        self.matrix[idx_snarl, idx_geno] = 1

    def __str__(self) :
        return f"{self.matrix[0]} \n" \
               f"{self.matrix[1]} \n" \
               f"{self.matrix[2]} \n" \
               f"{self.matrix[3]} \n" \
               f"{self.matrix[4]} \n"

class SnarlProcessor:
    def __init__(self, vcf_path: str):
        self.list_samples = VCF(vcf_path).samples
        self.matrix = Matrix(1_000_000, len(self.list_samples)*2)
        self.vcf_path = vcf_path

    def expand_matrix(self):
        """Expands a given numpy matrix by doubling the number of rows."""
        data_matrix = self.matrix.get_matrix()
        current_rows, current_cols = data_matrix.shape
        new_rows = current_rows + self.matrix.get_default_row_number() # add + default_row_number row  

        # Create a new matrix of zeros with the expanded size
        expanded_matrix = np.zeros((new_rows, current_cols), dtype=data_matrix.dtype)
        expanded_matrix[:current_rows, :] = data_matrix
        self.matrix.set_matrix(expanded_matrix)
            
    def determine_str(self, s: str, length_s : int, i: int) -> tuple[int, int]:
        """Extract an integer from a string starting at index i."""
        start_idx = i
        while i < length_s and s[i] not in ['>', '<']:
            i += 1
        return i, s[start_idx:i]

    def decompose_string(self, s: str) -> List[str]:
        """Decompose a string with snarl information."""
        result = []
        i = 0
        length_s = len(s)
        prev_int = None
        prev_sym = None
        
        while i < length_s:
            start_sym = s[i]
            i += 1
            i, current_int = self.determine_str(s, length_s, i)

            if prev_int is not None and prev_sym is not None:
                result.append(f"{prev_sym}{prev_int}{start_sym}{current_int}")
            
            prev_int = current_int
            prev_sym = start_sym
        
        return result
        
    def decompose_snarl(self, lst: List[str]) -> List[List[str]]:
        """Decompose a list of snarl strings."""
        return [self.decompose_string(s) for s in lst]

    def get_or_add_index(self, ordered_dict, key, length_ordered_dict):
        """ 
        Retrieve the index of the key if it exists in the OrderedDict.
        If the key does not exist, add it and return the new index.
        """
        if key in ordered_dict:
            return ordered_dict[key]
        else:
            new_index = length_ordered_dict
            ordered_dict[key] = new_index
            return new_index
    
    def push_matrix(self, idx_snarl, decomposed_snarl, row_header_dict, index_column):
        """Add True to the matrix if snarl is found"""

        # Retrieve or add the index in one step and calculate the length once
        length_ordered_dict = len(row_header_dict)
        idx_snarl = self.get_or_add_index(row_header_dict, decomposed_snarl, length_ordered_dict)

        # Check if a new matrix chunk is needed (only if length > 1)
        current_rows_number = self.matrix.get_matrix().shape[0]
        if length_ordered_dict > current_rows_number -1 :
            self.expand_matrix()

        # Add data to the matrix
        self.matrix.add_data(idx_snarl, index_column)

    def fill_matrix(self):
        """Parse VCF file (main function)"""
        row_header_dict = dict()

        # Parse variant line by line
        for variant in VCF(self.vcf_path):
            genotypes = variant.genotypes  # Extract genotypes once per variant
            snarl_list = variant.INFO.get('AT', '').split(',')  # Extract and split snarl list once per variant
            list_list_decomposed_snarl = self.decompose_snarl(snarl_list)  # Decompose snarls once per variant

            for index_column, genotype in enumerate(genotypes) :

                allele_1, allele_2 = genotype[:2]  # assume there are only 2 allele
                col_idx = index_column * 2

                if allele_1 == -1 or allele_2 == -1 : # case where we got ./.
                    continue

                for decompose_allele_1 in list_list_decomposed_snarl[allele_1] :
                    self.push_matrix(allele_1, decompose_allele_1, row_header_dict, col_idx)

                for decompose_allele_2 in list_list_decomposed_snarl[allele_2] :
                    self.push_matrix(allele_2, decompose_allele_2, row_header_dict, col_idx + 1)

        self.matrix.set_row_header(row_header_dict)

    def check_pheno_group(self, group) :
        """Check if all sample name in the matrix are matching with phenotype else return error"""
        if type(group) == tuple :
            list_group = list(group[0].keys()) + list(group[1].keys())
        elif type(group) == dict :
            list_group = [i for i in group.keys()]
        else :
            raise ValueError(f"group type : {type(group)} not an dict or a tuple.")

        set_sample = set(self.list_samples)
        set_group = set(list_group)
        missing_elements = set_sample - set_group

        if missing_elements:
            raise ValueError(f"The following sample name from merged vcf are not present in group file : {missing_elements}")

    def binary_table(self, snarls, binary_groups, output="output/binary_output.tsv") : 

        self.check_pheno_group(binary_groups)

        with open(output, 'wb') as outf:
            headers = 'CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tP_Fisher\tP_Chi2\tTable_sum\tNumber_column\tInter_group\tAverage\n'
            outf.write(headers.encode('utf-8'))

            for snarl, list_snarl in snarls.items() :
                df = self.create_binary_table(binary_groups, list_snarl)
                fisher_p_value, chi2_p_value, total_sum, numb_colum, inter_group, average = self.binary_stat_test(df)
                chrom = pos = type_var = ref = alt = ""
                data = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, snarl, type_var, ref, alt, fisher_p_value, chi2_p_value, total_sum, numb_colum, inter_group, average)
                outf.write(data.encode('utf-8'))
 
    def quantitative_table(self, snarls, quantitative, output="output/quantitative_output.tsv") :

        self.check_pheno_group(quantitative)

        with open(output, 'wb') as outf:
            headers = 'CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tP\tBETA\tSE\n'
            outf.write(headers.encode('utf-8'))
            for snarl, list_snarl in snarls.items() :
                df = self.create_quantitative_table(list_snarl)
                p_value, beta, se = self.linear_regression(df, quantitative)
                chrom = pos = type_var = ref = alt = ""
                data = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, snarl, type_var, ref, alt, p_value, beta, se)
                outf.write(data.encode('utf-8'))

    def identify_correct_path(self, decomposed_snarl: list, row_headers_dict: dict, idx_srr_save: list) -> list:
        """
        Return a list of column index where all specifique element of this column of matrix are 1
        """
        matrix = self.matrix.get_matrix()
        rows_to_check = np.array([], dtype=int)

        for snarl in decomposed_snarl :
            if "*" in snarl :
                continue
            if snarl in row_headers_dict :
                rows_to_check = np.append(rows_to_check, row_headers_dict[snarl])
            else:
                return [] 
        
        extracted_rows = matrix[rows_to_check, :]
        columns_all_ones = np.all(extracted_rows == 1, axis=0)
        idx_srr_save = np.where(columns_all_ones)[0].tolist()

        return idx_srr_save

    def create_binary_table(self, groups, list_path_snarl) -> pd.DataFrame :
        """Generates a binary table DataFrame indicating the presence of snarl paths in given groups based on matrix data"""
        row_headers_dict = self.matrix.get_row_header()
        list_samples = self.list_samples
        length_column_headers = len(list_path_snarl)

        # Initialize g0 and g1 with zeros, corresponding to the length of column_headers
        g0 = [0] * length_column_headers
        g1 = [0] * length_column_headers

        # Iterate over each path_snarl in column_headers
        for idx_g, path_snarl in enumerate(list_path_snarl):
            idx_srr_save = list(range(len(list_samples)))
            decomposed_snarl = self.decompose_string(path_snarl)
            idx_srr_save = self.identify_correct_path(decomposed_snarl, row_headers_dict, idx_srr_save)
            
            # Count occurrences in g0 and g1 based on the updated idx_srr_save
            for idx in idx_srr_save :
                srr = list_samples[idx//2]

                if srr in groups[0]:
                    g0[idx_g] += 1
                if srr in groups[1]:  
                    g1[idx_g] += 1

        # Create and return the DataFrame
        df = pd.DataFrame([g0, g1], index=['G0', 'G1'], columns=list_path_snarl)
        return df
    
    def create_quantitative_table(self, column_headers: list) -> pd.DataFrame:
        row_headers_dict = self.matrix.get_row_header()
        length_sample = len(self.list_samples)

        # Initialize a zero matrix for genotypes with shape (length_sample, len(column_headers))
        genotypes = np.zeros((length_sample, len(column_headers)), dtype=int)

        # Iterate over each path_snarl and fill in the matrix
        for col_idx, path_snarl in enumerate(column_headers):
            decomposed_snarl = self.decompose_string(path_snarl)
            idx_srr_save = self.identify_correct_path(decomposed_snarl, row_headers_dict, list(range(length_sample)))

            for idx in idx_srr_save:
                srr_idx = idx // 2  # Convert index to the appropriate sample index
                genotypes[srr_idx, col_idx] += 1

        df = pd.DataFrame(genotypes, index=self.list_samples, columns=column_headers)
        return df

    def sm_ols(self, x, y) :

        x_with_const = sm.add_constant(x)
        result = sm.OLS(y, x_with_const).fit()
        beta = result.params  # The coefficients (including intercept)
        se = result.bse  # The standard errors of the coefficients
        p_value = result.f_pvalue
        # print(result.summary())

        return p_value, beta, se
    
    def linear_regression(self, df, pheno : dict) -> float :
        
        df = df.astype(int)
        df['Target'] = df.index.map(pheno)
        x = df.drop('Target', axis=1)
        y = df['Target']
        p_value, beta, se = self.sm_ols(x, y)

        return p_value, beta, se

    def chi2_test(self, df) -> float:
        """Calculate p_value from list of dataframe using chi-2 test"""

        try:
            # Perform Chi-Square test
            p_value = chi2_contingency(df)[1]

        except ValueError as e:
            p_value = "Error"
    
        return p_value

    def fisher_test(self, df) -> float : 
        """Calcul p_value using fisher exact test"""

        try:
            # Perform Fisher exact test
            p_value = fisher_exact(df)[1]

        except ValueError as e: 
            p_value = "Error"
        
        return p_value
     
    def binary_stat_test(self, df) :
 
        fisher_p_value = self.fisher_test(df)
        chi2_p_value = self.chi2_test(df)
        total_sum = int(df.values.sum())
        inter_group = int(df.min().sum())
        numb_colum = df.shape[1]
        average = float(total_sum / numb_colum)
        return fisher_p_value, chi2_p_value, total_sum, numb_colum, inter_group, average

def parse_group_file(group_file : str):

    df = pd.read_csv(group_file, sep='\t')
    
    # Create dictionaries for group 0 and group 1
    group_0 = {sample: 0 for sample in df[df['GROUP'] == 0]['SAMPLE']}
    group_1 = {sample: 1 for sample in df[df['GROUP'] == 1]['SAMPLE']}
    return group_0, group_1
 
def parse_pheno_file(file_path : str) -> dict:

    df = pd.read_csv(file_path, sep='\t')

    # Extract the IID (second column) and PHENO (third column) and convert PHENO to float
    parsed_pheno = dict(zip(df['IID'], df['PHENO']))
    return parsed_pheno

def parse_snarl_path_file(path_file: str) -> dict:
    
    # Initialize an empty dictionary for the snarl paths
    snarl_paths = defaultdict(list)

    # Read the file into a pandas DataFrame
    df = pd.read_csv(path_file, sep='\t', dtype=str)
    df = df[df['paths'].notna()]
    df['paths'] = df['paths'].str.split(',')

    # Create the dictionary with snarl as key and paths as values
    for snarl, paths in zip(df['snarl'], df['paths']):
        snarl_paths[snarl].extend(paths)

    return snarl_paths

def check_format_vcf_file(file_path : str) -> str:
    """
    Function to check if the provided file path is a valid VCF file.
    """
    if not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")

    if not file_path.lower().endswith('.vcf') and not file_path.lower().endswith('.vcf.gz'):
        raise argparse.ArgumentTypeError(f"The file {file_path} is not a valid VCF file. It must have a .vcf extension or .vcf.gz.")
    return file_path

def check_format_group_snarl(file_path : str) -> str :
    """
    Function to check if the provided file path is a valid group file.
    """
    if not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")
    
    if not file_path.lower().endswith('.txt') and not file_path.lower().endswith('.tsv'):
        raise argparse.ArgumentTypeError(f"The file {file_path} is not a valid group/snarl file. It must have a .txt extension or .tsv.")
    return file_path
    
def check_format_pheno_q(file_path : str) -> str :

    if not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")
    
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
    
    header = first_line.split('\t')
    expected_header = ['FID', 'IID', 'PHENO']
    if header != expected_header:
        raise ValueError(f"The file must contain the following headers: {expected_header} and be split by tabulation")
    
    return file_path

def check_format_pheno_b(file_path : str) -> str :

    if not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")
    
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
    
    header = first_line.split('\t')
    expected_header = ['SAMPLE', 'GROUP']
    if header != expected_header:
        raise ValueError(f"The file must contain the following headers: {expected_header} and be split by tabulation")
    
    return file_path

def get_first_snarl(s):

    match = re.findall(r'\d+', s)
    if match:
        return int(match[0])
    return None  # Return None if no integers are found

def classify_variant(ref, alt) :

    if len(ref) == len(alt) == 1:
        return "SNP"
    elif len(ref) > len(alt) and alt != "<DEL>":
        return "DEL"
    elif len(ref) < len(alt):
        return "INS"
    elif len(ref) == len(alt) and len(ref) > 1:
        return "MNP"
    else :
        raise ValueError(f"what is this ref : {ref}, alt : {alt}")

def write_pos_snarl(vcf_file, output_file):
    vcf_dict = parse_vcf_to_dict(vcf_file)
    
    # Read the output file, fill placeholders, and collect lines for rewrite
    with open(output_file, 'r', encoding='utf-8') as out_f:
        lines = out_f.readlines()

    with open(output_file, 'w', encoding='utf-8') as out_f:
        for line in lines:
            columns = line.strip().split('\t')
            snarl = columns[2]  # Assuming SNARL is in column 3 (index 2)
            info = match_pos(snarl, vcf_dict)
            
            if info:
                chrom, pos, type_var, ref, alt = info
            else:
                chrom = pos = type_var = ref = alt = "NA"
            
            # Replace placeholders with actual values in the correct columns
            columns[0] = chrom
            columns[1] = pos
            columns[3] = type_var
            columns[4] = ref
            columns[5] = alt

            # Write the modified line
            out_f.write('\t'.join(columns) + '\n')

def match_pos(snarl, vcf_dict):
    """Matches the SNARL to an entry in the VCF dictionary, if available."""
    start_snarl, _ = snarl.split('_')
    return vcf_dict.get(start_snarl, None)

def parse_vcf_to_dict(vcf_file):
    """Parses a VCF file and returns a dictionary with SNARL IDs as keys."""
    vcf_dict = {}
    
    for record in VCF(vcf_file):
        chrom = record.CHROM  # Chromosome
        pos = record.POS      # Position
        snarl = get_first_snarl(record.ID) if record.ID else None
        ref = record.REF      # Reference allele
        alt = record.ALT[0]   # First alternative allele (assuming biallelic)
        
        variant_type = classify_variant(ref, alt)
        if snarl:  # Only add entries with a valid snarl
            vcf_dict[snarl] = (chrom, pos, variant_type, ref, alt)

    return vcf_dict

if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description="Parse and analyse snarl from vcf file")
    parser.add_argument("vcf_path", type=check_format_vcf_file, help="Path to the vcf file (.vcf or .vcf.gz)")
    parser.add_argument("snarl", type=check_format_group_snarl, help="Path to the snarl file that containt snarl and aT (.txt or .tsv)")
    parser.add_argument("vcf_pangenome", type=check_format_vcf_file, help="Path to the vcf referencing all position of snarl")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--binary", type=check_format_pheno_b, help="Path to the binary group file (.txt or .tsv)")
    group.add_argument("-q", "--quantitative", type=check_format_pheno_q, help="Path to the quantitative phenotype file (.txt or .tsv)")
    parser.add_argument("-o", "--output", type=str, required=False, help="Path to the output file")
    args = parser.parse_args()
    
    start = time.time()
    vcf_object = SnarlProcessor(args.vcf_path)
    vcf_object.fill_matrix()
    print(f"Time Matrix : {time.time() - start} s")

    start = time.time()
    snarl = parse_snarl_path_file(args.snarl)
    output = args.output if args.output else None

    if args.binary:
        binary_group = parse_group_file(args.binary)
        vcf_object.binary_table(snarl, binary_group, output)

    if args.quantitative:
        quantitative = parse_pheno_file(args.quantitative)
        vcf_object.quantitative_table(snarl, quantitative, output)

    #write_pos_snarl(args.vcf_pangenome, args.output)
    print(f"Time : {time.time() - start} s")

