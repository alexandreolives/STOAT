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
import time

class Chunk:
    def __init__(self, MATRIX_ROW_NUMBER, nb_matrix_column):
        self.data = np.zeros((MATRIX_ROW_NUMBER, nb_matrix_column), dtype=bool)

    def add_data(self, idx_snarl, idx_geno):
        self.data[idx_snarl, idx_geno] = 1

    def get_data(self) :
        return self.data
    
    def __str__(self) :
        return f"{self.data}"
    
    @staticmethod
    def concatenate_matrix(matrix_list, axis=0):
        
        matrix_list = [matrix.get_data() for matrix in matrix_list]
        if len(matrix_list) == 0 :
            raise ValueError("The array list is empty, nothing to concatenate.")
        
        if len(matrix_list) == 1 :
            return matrix_list[0]
        
        # Check if all matrix have the same shape except for the concatenation axis
        reference_shape = list(matrix_list[0].shape)
        reference_shape[axis] = None  # Ignore the axis we're concatenating along
        for matrix in matrix_list:
            if list(matrix.shape)[:axis] + list(matrix.shape)[axis+1:] != reference_shape[:axis] + reference_shape[axis+1:]:
                raise ValueError("All arrays must have the same shape except along the concatenation axis.")

        return np.concatenate(matrix_list, axis=0)

class Matrix :
    def __init__(self, Tupple_matrix : tuple) :
        self.matrix = Tupple_matrix[0]
        self.row_header = Tupple_matrix[1]
        self.column_header = Tupple_matrix[2]
    
    def get_matrix(self) :
        return self.matrix
    
    def get_row_header(self) :
        return self.row_header
    
    def get_column_header(self) :
        return self.column_header
    
    def __str__(self) :
        return f"           {self.column_header} \n" \
               f"{self.row_header[0]} {self.matrix[0]}"
        
class Snarl :
    def __init__(self, vcf_path: str) :
        self.vcf_path = vcf_path
        self.matrix = None
        self.list_samples = []
    
    def initialise_matrix(self) :
        self.matrix = Matrix(self.create_matrix())

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

    def push_matrix(self, decomposed_snarl, allele, MATRIX_ROW_NUMBER, nb_matrix_column, row_header, row_header_dict, Matrix_list, idx_geno, idx_matrix):
        """Add True to the matrix if snarl are found"""

        current_row_len = len(row_header)

        for decomposed_snarl_path in decomposed_snarl[allele]:
            # Check if the path is already in the row header dictionary
            if decomposed_snarl_path not in row_header_dict:
                row_header_dict[decomposed_snarl_path] = current_row_len
                row_header.append(decomposed_snarl_path)
                current_row_len += 1 

            idx_snarl = row_header_dict[decomposed_snarl_path]

            # Determine if a new matrix chunk is needed
            if current_row_len % MATRIX_ROW_NUMBER == 1:
                # Add new chunk
                Matrix_list.append(Chunk(MATRIX_ROW_NUMBER, nb_matrix_column))
                idx_matrix += 1

            # Use the current matrix to add data
            Matrix_list[idx_matrix].add_data(idx_snarl % MATRIX_ROW_NUMBER, idx_geno) # matrix[idx_snarl, idx_geno] = 1

        return Matrix_list, idx_matrix

    def create_matrix(self):
        """Parse VCF file (main function)"""
        vcf_object = VCF(self.vcf_path)
        self.list_samples = vcf_object.samples
        # define column and row header for the numpy matrix 
        row_header = []
        row_header_dict = {}  # Use dictionary for faster lookup
        column_header = [f"{sample}_{i}" for sample in self.list_samples for i in range(2)]
        nb_matrix_column = len(self.list_samples) * 2
        MATRIX_ROW_NUMBER = 100000
        Matrix_list = [Chunk(MATRIX_ROW_NUMBER, nb_matrix_column)]
        idx_matrix = 0

        # Parse variant line per line
        for idx, variant in enumerate(vcf_object):
            genotypes = variant.genotypes  # genotypes: [[-1, -1, False], [-1, -1, False], [1, 1, False], [-1, -1, False]]
            snarl_list = variant.INFO.get('AT', '').split(',')  # snarl_list: ['>1272>1273>1274', '>1272>1274']
            decomposed_snarl = self.decompose_snarl(snarl_list)  # decomposed_snarl: [['>1272>1273', '>1273>1274'], ['>1272>1274']]
            idx_geno = 0

            for gt in genotypes:
                allele_1, allele_2 = gt[:2]  # Extract the two alleles

                if allele_1 != -1 and allele_2 != -1:  # TODO SPECIAL CASE 
                    Matrix_list, idx_matrix = self.push_matrix(decomposed_snarl, allele_1, MATRIX_ROW_NUMBER, nb_matrix_column, row_header, row_header_dict, Matrix_list, idx_geno, idx_matrix)
                    Matrix_list, idx_matrix = self.push_matrix(decomposed_snarl, allele_2, MATRIX_ROW_NUMBER, nb_matrix_column, row_header, row_header_dict, Matrix_list, idx_geno + 1, idx_matrix)

                idx_geno += 2  # push the index

        print(f"Decomposed snarl found : {len(row_header)}")
        return Chunk.concatenate_matrix(Matrix_list), row_header, column_header, self.list_samples
    
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
        for reference, snarl in snarls.items() :
            df = self.create_binary_table(binary_groups, snarl)
            res_table.append(df)
            reference_list.append(reference)

        return reference_list, res_table
    
    def quantitative_table(self, snarls, pheno) -> list : # list of dataframe
        self.check_pheno_group(pheno)
        res_table = []
        res_reference = []
        for ref, snarl in snarls.items() :
            df = self.create_quantitative_table(snarl)
            res_table.append(df)
            res_reference.append(ref)

        return res_reference, res_table
    
    def create_quantitative_table(self, column_headers) -> pd.DataFrame:

        #print("column_headers : ", column_headers)
        row_headers = self.matrix.get_row_header()
        row_headers_dict = {header: idx for idx, header in enumerate(row_headers)}  # Dictionary for faster lookup
        column_headers_header = self.list_samples
        matrix = self.matrix.get_matrix()  # Retrieve matrix once
        genotypes = []

        # Iterate over each path_snarl in column_headers
        for path_snarl in column_headers:
            idx_srr_save = set(range(len(column_headers_header)))  # Use a set for efficient filtering
            decomposed_snarl = self.decompose_string(path_snarl)

            for snarl in decomposed_snarl:
                if "*" in snarl:
                    continue

                if snarl in row_headers_dict:
                    idx_row = row_headers_dict[snarl]
                    matrix_row = matrix[idx_row]

                    # Efficiently update idx_srr_save with indices where matrix_row is 1
                    idx_srr_save.intersection_update(idx for idx in idx_srr_save if any(matrix_row[2 * idx: 2 * idx + 2]))

                    # Early exit if no valid indices are left
                    if not idx_srr_save:
                        break
                else:
                    idx_srr_save.clear()  # Clear the set, indicating no match found
                    break

            genotype = [1 if idx in idx_srr_save else 0 for idx in range(len(column_headers_header))]
            genotypes.append(genotype)

        # Transposing the matrix
        transposed_genotypes = list(map(list, zip(*genotypes)))
        print("transposed_genotypes : ", transposed_genotypes)
        print("column_headers_header : ", column_headers_header)

        df = pd.DataFrame(transposed_genotypes, index=column_headers_header, columns=column_headers)
        
        return df

    def create_binary_table(self, groups, column_headers) -> pd.DataFrame :
        
        row_headers = self.matrix.get_row_header()
        column_headers_header = self.list_samples
        length_column_headers = len(column_headers)

        # Initialize g0 and g1 with zeros, corresponding to the length of column_headers
        g0 = [0] * length_column_headers
        g1 = [0] * length_column_headers

        # Iterate over each path_snarl in column_headers
        for idx_g, path_snarl in enumerate(column_headers):
            idx_srr_save = list(range(len(column_headers_header)))  # Initialize with all indices
            decomposed_snarl = self.decompose_string(path_snarl)
            # Iterate over each snarl in the decomposed_snarl
            for snarl in decomposed_snarl:
                
                # Suppose that at leat one snarl pass thought
                # Case * in snarl
                if "*" in snarl:
                    continue

                else :
                    if any(snarl in header for header in row_headers) :
                        idx_row = row_headers.index(snarl)
                        matrix_row = self.matrix.get_matrix()[idx_row]

                        # Update idx_srr_save only where indices row is 1
                        idx_srr_save = [idx for idx in idx_srr_save if any(matrix_row[2 * idx: 2 * idx + 2])]
                                        
                    else :
                        idx_srr_save = []
                        break

            # Count occurrences in g0 and g1 based on the updated idx_srr_save
            for idx in idx_srr_save:
                srr = column_headers_header[idx]
                if srr in groups[0]:
                    g0[idx_g] += 1
                if srr in groups[1]:  
                    g1[idx_g] += 1

        # Create and return the DataFrame
        df = pd.DataFrame([g0, g1], index=['G0', 'G1'], columns=column_headers)
        #print(df)
        return df

    def linear_regression(self, list_dataframe, pheno) :

        list_pvalue = []
        for df in list_dataframe :
                
            df = df.astype(int)
            df['Target'] = df.index.map(pheno)

            x = df.drop('Target', axis=1)
            y = df['Target']

            # Fit the regression model
            model = sm.OLS(y, x).fit()

            # Extract p-values from the fitted model
            p_values = model.pvalues
            #print("p_values : ", p_values)
            list_pvalue.append(p_values)
        return list_pvalue

    def chi2_test(self, list_dataframe) -> list:
        """Calculate p_value from list of dataframe using chi-2 test"""

        list_pvalue = []
        for df in list_dataframe:
            #print("df : ", df)
            
            # Check if dataframe has at least 2 columns and more than 0 counts in every cell
            if df.shape[1] >= 2 and np.all(df.sum(axis=0)) and np.all(df.sum(axis=1)):
                try:
                    # Perform Chi-Square test
                    chi2, p_value, dof, expected = chi2_contingency(df)
                    #print(f"Chi-Square Test p-value: {p_value}")
                except ValueError as e:
                    print(f"Error in Chi-Square test: {e}")
                    p_value = "Error"
            else:
                p_value = "N/A"

            list_pvalue.append(p_value)
            
        return list_pvalue

    def fisher_test(self, list_dataframe) : 
        """Calcul p_value from list of dataframe using fisher exact test"""

        list_pvalue = []
        for df in list_dataframe :
            #print(df)
            try:
                odds_ratio, p_value = fisher_exact(df)
            except ValueError as e: 
                p_value = 'N/A'

            list_pvalue.append(p_value)
            #print(f"Fisher's Exact Test p-value: {p_value}")
        
        return list_pvalue
     
    def binary_stat_test(self, list_dataframe) :

        fisher_p_value = self.fisher_test(list_dataframe)
        chi2_p_value = self.chi2_test(list_dataframe)

        return fisher_p_value, chi2_p_value

    def output_writing_binary(self, reference_list, list_pvalues, output_filename="binary_output.tsv") :
        # Combine DataFrames and p-values into a single DataFrame for saving
        list_ficher = list_pvalues[0]
        list_chi = list_pvalues[1]
        df_combined = pd.DataFrame({
            'Snarl': reference_list,
            'P_value (Fisher)': list_ficher,
            'P_value (Chi2)': list_chi
        })
        
        # Save to TSV
        df_combined.to_csv(output_filename, sep='\t', index=False)

    def output_writing_quantitative(self, reference_list, list_pvalues, output_filename="quantitative_output.tsv") :
        print("reference_list : ", reference_list)
        df_combined = pd.DataFrame({
                'Snarl': reference_list,
                'P_value': list_pvalues
            })

        # Save to TSV
        df_combined.to_csv(output_filename, sep='\t', index=False)

def parse_group_file(groupe_file):
    group_0 = []
    group_1 = []

    with open(groupe_file, 'r') as file:
        for line in file:
            if line.strip():  # Check if the line is not empty
                parts = line.split()
                group_0.append(parts[0])

                if len(parts) > 1:
                    group_1.append(parts[1])
                else:
                    group_1.append('')  # Append an empty string if the second group is missing

    return group_0, group_1

def parse_pheno_file(file_path):

    parsed_data = {}
    with open(file_path, 'r') as file:
        # Skip the header line
        next(file)
        
        # Process each subsequent line
        for line in file:
            parts = line.strip().split()
            if len(parts) != 3:
                raise ValueError(f"The line '{line.strip()}' does not have exactly 3 columns.")
            
            sample_name = parts[1]  # IID
            phenotype = parts[2]    # PHENO
            parsed_data[sample_name] = float(phenotype)
    
    return parsed_data

def parse_snarl_path_file(path_file: str) -> dict:
    path_list = defaultdict(list)

    with open(path_file, 'r') as file:
        next(file)

        for line in file:
            if line.strip():  # Skip empty lines
                snarl, aT = line.split()  
                path_list[snarl].append(aT)

    return path_list

def check_format_vcf_file(file_path):
    """
    Function to check if the provided file path is a valid VCF file.
    """
    if not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")

    if not file_path.lower().endswith('.vcf') and not file_path.lower().endswith('.vcf.gz'):
        raise argparse.ArgumentTypeError(f"The file {file_path} is not a valid VCF file. It must have a .vcf extension or .vcf.gz.")
    return file_path

def check_format_group_snarl(file_path):
    """
    Function to check if the provided file path is a valid group file.
    """
    if not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")
    
    if not file_path.lower().endswith('.txt') and not file_path.lower().endswith('.tsv'):
        raise argparse.ArgumentTypeError(f"The file {file_path} is not a valid group/snarl file. It must have a .txt extension or .tsv.")
    return file_path
    
def check_format_pheno(file_path):

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
    vcf_object = Snarl(args.vcf_path)
    vcf_object.initialise_matrix()
    print(f"Time : {time.time() - start}")
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
        reference_list, list_quantitative_df = vcf_object.quantitative_table(snarl, quantitative)
        quantitative_p_value = vcf_object.linear_regression(list_quantitative_df, quantitative)
        if args.output :
            vcf_object.output_writing_quantitative(reference_list, quantitative_p_value, args.output)
        else :
            vcf_object.output_writing_quantitative(reference_list, quantitative_p_value)

    print(f"Time : {time.time() - start} s")

