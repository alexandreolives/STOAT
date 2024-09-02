import argparse
from cyvcf2 import VCF
from typing import List
import numpy as np
import pandas as pd
import statsmodels.api as sm
from collections import defaultdict
from collections import OrderedDict
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
import os
import time
    
class Matrix :
    def __init__(self, default_row_number=100000, column_number=2):
        self.matrix = np.zeros((default_row_number, column_number),dtype=bool)
        self.row_header = OrderedDict()
    
    def get_matrix(self):
        return self.matrix

    def get_row_header(self):
        return self.row_header

    def set_row_header(self, row_header):
        self.row_header = row_header  

    def add_data(self, idx_snarl, idx_geno):
        self.matrix[idx_snarl, idx_geno] = 1

    def __str__(self) :
        return f"{self.row_header[0]} {self.matrix[0]} \n" \
               f"{self.row_header[1]} {self.matrix[1]} \n" \
               f"{self.row_header[2]} {self.matrix[2]} \n" \
               f"{self.row_header[3]} {self.matrix[3]} \n" \
               f"{self.row_header[4]} {self.matrix[4]} \n"

class SnarlProcessor:
    def __init__(self, vcf_path: str):
        self.list_samples = VCF(vcf_path).samples
        self.matrix = Matrix(100000, len(self.list_samples)*2)
        self.vcf_path = vcf_path

    def expand_matrix(self):
        """
        Expands a given numpy matrix by doubling the number of rows.
        """

        data_matrix = self.matrix.get_matrix()
        current_rows, current_cols = data_matrix.shape
        new_rows = current_rows * 2

        # Create a new matrix of zeros with the expanded size
        expanded_matrix = np.zeros((new_rows, current_cols), dtype=data_matrix.dtype)
        expanded_matrix[:current_rows, :] = data_matrix

        return expanded_matrix
            
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
    
    def push_matrix(self, idx_snarl, decomposed_snarl, allele, row_header_dict, index_column):
        """Add True to the matrix if snarl is found"""

        current_rows, _ = self.matrix.get_matrix().shape

        if allele == idx_snarl:
            # Retrieve or add the index in one step and calculate the length once
            length_ordered_dict = len(row_header_dict)
            idx_snarl = self.get_or_add_index(row_header_dict, decomposed_snarl, length_ordered_dict)

            # Check if a new matrix chunk is needed (only if length > 1)
            if length_ordered_dict > 1 and length_ordered_dict > current_rows :
                self.expand_matrix()

            # Compute index row
            row_index = idx_snarl
            
            # Add data to the matrix
            self.matrix.add_data(row_index, index_column)

    def truncate_matrix(self, max_rows):
        """
        Truncates the matrix to a specified maximum number of rows.
        """

        data_matrix = self.matrix.get_matrix()
        current_rows, _ = data_matrix.shape
        if max_rows == current_rows:
            return data_matrix

        truncated_matrix = data_matrix[:max_rows, :]

        return truncated_matrix

    def fill_matrix(self): # TODO : special case where snarl decomposed are between 2 chunk ... 
        """Parse VCF file (main function)"""
        row_header_dict = OrderedDict()

        # Parse variant line by line
        for variant in VCF(self.vcf_path):
            genotypes = variant.genotypes  # Extract genotypes once per variant
            snarl_list = variant.INFO.get('AT', '').split(',')  # Extract and split snarl list once per variant
            list_list_decomposed_snarl = self.decompose_snarl(snarl_list)  # Decompose snarls once per variant

            # Loop over each decomposed snarl list and genotype
            for idx_snarl, list_decomposed_snarl in enumerate(list_list_decomposed_snarl):
                for decomposed_snarl in list_decomposed_snarl:
                    # Loop over genotypes with index_column tracking
                    for index_column, gt in enumerate(genotypes):
                        allele_1, allele_2 = gt[:2]  # Extract alleles
                        
                        # Push matrix only if both alleles are valid
                        if allele_1 != -1 and allele_2 != -1:
                            col_idx = index_column * 2
                            self.push_matrix(idx_snarl, decomposed_snarl, allele_1, row_header_dict, col_idx)
                            self.push_matrix(idx_snarl, decomposed_snarl, allele_2, row_header_dict, col_idx + 1)

        self.matrix.set_row_header(row_header_dict)

    def check_pheno_group(self, group) :
        """Check if all sample name in the matrix are matching with phenotype else return error"""
        if type(group) == tuple :
            list_group = group[0] + group[1]
        elif type(group) == dict :
            list_group = [i for i in group.keys()]
        else :
            raise ValueError(f"group type : {type(group)} not an dict or a tuple.")

        set_sample = set(self.list_samples)
        set_group = set(list_group)
        missing_elements = set_sample - set_group

        if missing_elements:
            raise ValueError(f"The following elements from set_sample are not present in set_group: {missing_elements}")

    def binary_table(self, snarls, binary_groups) -> list : # list of dataframe
        self.check_pheno_group(binary_groups)
        res_table = []
        reference_list = []
        for reference, list_snarl in snarls.items() :
            df = self.create_binary_table(binary_groups, list_snarl)
            res_table.append(df)
            reference_list.append(reference)

        return reference_list, res_table

    def identify_correct_path(self, decomposed_snarl: list, row_headers_dict: dict, idx_srr_save: list) -> list:
        """
        Return a list of column index where all specifique element of this column of matrix are 1
        """
        matrix = self.matrix.get_matrix()
        rows_to_check = []
        
        for snarl in decomposed_snarl :
            if "*" not in snarl and snarl in row_headers_dict :
                idx_row = row_headers_dict[snarl]
                rows_to_check.append(idx_row)
            else:
                return [] 
            
        if not rows_to_check :
            return []

        rows_to_check = np.array(rows_to_check)
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

    def quantitative_table(self, snarls, pheno) -> list : # list of dataframe
        self.check_pheno_group(pheno)
        res_table = []
        for _, snarl in snarls.items() :
            res_table.append(self.create_quantitative_table(snarl))

        return res_table
    
    def create_quantitative_table(self, column_headers : list) -> pd.DataFrame:
        row_headers_dict = self.matrix.get_row_header()
        column_headers_header = self.list_samples
        genotypes = []

        # Iterate over each path_snarl in column_headers
        for path_snarl in column_headers:
            idx_srr_save = list(range(len(column_headers_header)))
            decomposed_snarl = self.decompose_string(path_snarl)
            idx_srr_save = self.identify_correct_path(decomposed_snarl, row_headers_dict, idx_srr_save)

            genotype = [1 if idx in idx_srr_save else 0 for idx in range(len(column_headers_header))]
            genotypes.append(genotype)

        # Transposing the matrix
        transposed_genotypes = list(map(list, zip(*genotypes)))
        df = pd.DataFrame(transposed_genotypes, index=column_headers_header, columns=column_headers)
        return df

    def linear_regression(self, list_dataframe : list, pheno : dict) -> list :

        list_pvalue = []
        for df in list_dataframe :
                
            df = df.astype(int)
            df['Target'] = df.index.map(pheno)

            x = df.drop('Target', axis=1)
            y = df['Target']

            # Fit the regression model
            model = sm.OLS(y, x).fit()

            # Extract p-values from the fitted model and format as a list of tuples
            for index, pval in model.pvalues.items() :
                if str(pval) == 'nan' :
                    pval = "N/A"
                list_pvalue.append((index, pval))
        return list_pvalue

    def chi2_test(self, dataframe) -> list:
        """Calculate p_value from list of dataframe using chi-2 test"""

        list_pvalue = []
        for df in dataframe:
            
            # Check if dataframe has at least 2 columns and more than 0 counts in every cell
            if df.shape[1] >= 2 and np.all(df.sum(axis=0)) and np.all(df.sum(axis=1)):
                try:
                    # Perform Chi-Square test
                    chi2, p_value, dof, expected = chi2_contingency(df)
                except ValueError as e:
                    p_value = "Error"
            else:
                p_value = "N/A"

            list_pvalue.append(p_value)
            
        return list_pvalue

    def fisher_test(self, list_dataframe : list) -> list : 
        """Calcul p_value from list of dataframe using fisher exact test"""

        list_pvalue = []
        for df in list_dataframe :
            try:
                odds_ratio, p_value = fisher_exact(df)

            except ValueError as e: 
                p_value = 'N/A'

            list_pvalue.append(p_value)
        
        return list_pvalue
     
    def binary_stat_test(self, list_dataframe : list) :

        fisher_p_value = self.fisher_test(list_dataframe)
        chi2_p_value = self.chi2_test(list_dataframe)

        return fisher_p_value, chi2_p_value

    def output_writing_binary(self, reference_list, list_pvalues, output_filename="output/binary_output.tsv") :
        # Combine DataFrames and p-values into a single DataFrame for saving
        list_ficher = list_pvalues[0]
        list_chi = list_pvalues[1]
        df_combined = pd.DataFrame({
            'Snarl': reference_list,
            'P_value (Fisher)': list_ficher,
            'P_value (Chi2)': list_chi
            #add number path
            #sum de toute la table
            #sum g0 et sum g1 et min (G0 ^ G1)
            #sum // column
        })
        
        df_combined.to_csv(output_filename, sep='\t', index=False)

    def output_writing_quantitative(self, list_pvalues, output_filename="output/quantitative_output.tsv") :

        df_combined = pd.DataFrame(list_pvalues, columns=['Snarl', 'P_value'])
        df_combined.to_csv(output_filename, sep='\t', index=False)

def parse_group_file(group_file : str):
    # Read the file into a DataFrame
    df = pd.read_csv(group_file, sep='\t')
    
    # Check if the required columns are present
    required_headers = {'sample', 'group'}
    if not required_headers.issubset(df.columns):
        raise ValueError(f"Missing required columns. The file must contain the following columns: {required_headers}")
    
    # Extract samples belonging to group 0 and group 1
    group_0 = df[df['group'] == 0]['sample'].tolist()
    group_1 = df[df['group'] == 1]['sample'].tolist()
    
    return group_0, group_1

def parse_pheno_file(file_path : str) -> dict:
    # Read the file into a DataFrame
    df = pd.read_csv(file_path, sep='\t')

    # Check if the required headers are present
    required_headers = {'FID', 'IID', 'PHENO'}
    if not required_headers.issubset(df.columns):
        raise ValueError(f"The file must contain the following headers: {required_headers}")

    # Extract the IID (second column) and PHENO (third column) and convert PHENO to float
    parsed_pheno = dict(zip(df['IID'], df['PHENO']))

    return parsed_pheno

def parse_snarl_path_file(path_file: str) -> dict:
    # Initialize a defaultdict with lists as default values
    snarl_paths = defaultdict(list)
    
    # Read the file using pandas
    df = pd.read_csv(path_file, sep='\t', dtype=str)
    
    # Iterate through the DataFrame rows
    for _, row in df.iterrows():
        snarl = row['snarl']
        paths = row['paths']
        
        # Split paths by comma and add them to the defaultdict
        if pd.notna(paths):
            path_list = paths.split(',')
            snarl_paths[snarl].extend(path_list)
    
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
    
def check_format_pheno(file_path : str) -> str :

    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
    
    header = first_line.split()
    if header != ['FID', 'IID', 'PHENO']:
        raise ValueError("The file does not contain the correct header.")
    return file_path

if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description="Parse and analyse snarl from vcf file")
    parser.add_argument("vcf_path", type=check_format_vcf_file, help="Path to the vcf file (.vcf or .vcf.gz)")
    parser.add_argument("snarl", type=check_format_group_snarl, help="Path to the snarl file that containt snarl and aT (.txt or .tsv)")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--binary", type=check_format_group_snarl, help="Path to the binary group file (.txt or .tsv)")
    group.add_argument("-q", "--quantitative", type=check_format_pheno, help="Path to the quantitative phenotype file (.txt or .tsv)")
    parser.add_argument("-o", "--output", type=str, required=False, help="Path to the output file")
    args = parser.parse_args()
    
    start = time.time()
    vcf_object = SnarlProcessor(args.vcf_path)
    vcf_object.fill_matrix()
    print(f"Time Matrix : {time.time() - start} s")
    start = time.time()

    snarl = parse_snarl_path_file(args.snarl)

    if args.binary:
        binary_group = parse_group_file(args.binary)
        reference_list, list_binary_df = vcf_object.binary_table(snarl, binary_group)
        binary_p_value = vcf_object.binary_stat_test(list_binary_df)

        if args.output :
            vcf_object.output_writing_binary(reference_list, binary_p_value, args.output)
        else :
            vcf_object.output_writing_binary(reference_list, binary_p_value)
        
    if args.quantitative:
        quantitative = parse_pheno_file(args.quantitative)
        list_quantitative_df = vcf_object.quantitative_table(snarl, quantitative)
        quantitative_p_value = vcf_object.linear_regression(list_quantitative_df, quantitative)

        if args.output :
            vcf_object.output_writing_quantitative(quantitative_p_value, args.output)
        else :
            vcf_object.output_writing_quantitative(quantitative_p_value)

    print(f"Time P-value: {time.time() - start} s")

