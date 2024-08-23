import argparse
from cyvcf2 import VCF
from typing import List
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import mean_squared_error, r2_score
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
    
    def decompose_snarl(self, lst: List[str]) -> List[List[str]] :
        """Decompose a list of snarl strings."""
        return [self.decompose_string(s) for s in lst]

    def push_matrix(self, decomposed_snarl, allele, MATRIX_ROW_NUMBER, nb_matrix_column, row_header, Matrix_list, idx_geno, idx_matrix) :
        """Add True to the matrix if snarl are found"""
        for decomposed_snarl_path in decomposed_snarl[allele] :
            if decomposed_snarl_path not in row_header :
                row_header.append(decomposed_snarl_path)
            
            idx_snarl = row_header.index(decomposed_snarl_path)
            
            # Check if current numpy matrix are full 
            if (len(row_header)+1) % MATRIX_ROW_NUMBER == 1 :
                # Add new chunck
                Matrix_list.append(Chunk(MATRIX_ROW_NUMBER, nb_matrix_column))
                idx_matrix += 1

            Matrix_list[idx_matrix].add_data(idx_snarl%MATRIX_ROW_NUMBER, idx_geno) # matrix[idx_snarl, idx_geno] = 1
        
        return Matrix_list, idx_matrix

    def create_matrix(self) :
        """Parse vcf file (main function)"""
        vcf_object = VCF(self.vcf_path)
        self.list_samples = vcf_object.samples
        # define column and row header for the numpy matrix 
        row_header = []
        column_header = [f"{sample}_{i}" for sample in self.list_samples for i in range(2)]
        nb_matrix_column = len(self.list_samples) * 2
        MATRIX_ROW_NUMBER = 100000
        Matrix_list = [Chunk(MATRIX_ROW_NUMBER, nb_matrix_column)]
        idx_matrix = 0

        # Parse variant line per line
        for variant in vcf_object :
            genotypes = variant.genotypes # genotypes : [[-1, -1, False], [-1, -1, False], [1, 1, False], [-1, -1, False]]
            snarl_list = variant.INFO.get('AT', '').split(',') # snarl_list :  ['>1272>1273>1274', '>1272>1274']
            decomposed_snarl = self.decompose_snarl(snarl_list) # decomposed_snarl :  [['>1272>1273', '>1273>1274'], ['>1272>1274']]

            idx_geno = 0
            for gt in genotypes:
                allele_1, allele_2 = gt[:2]  # Extract the two allele

                if allele_1 != -1 and allele_2 != -1 :  #TODO SPECIAL CASE 
                    Matrix_list, idx_matrix = self.push_matrix(decomposed_snarl, allele_1, MATRIX_ROW_NUMBER, nb_matrix_column, row_header, Matrix_list, idx_geno, idx_matrix)
                    Matrix_list, idx_matrix = self.push_matrix(decomposed_snarl, allele_2, MATRIX_ROW_NUMBER, nb_matrix_column, row_header, Matrix_list, idx_geno+1, idx_matrix)

                idx_geno += 2 # push the index

        return Chunk.concatenate_matrix(Matrix_list), row_header, column_header, self.list_samples

    def check_pheno_group(self, group) :
        """Check if all sample name in the matrix are matching with phenotype else return error"""
        if type(group) == tuple :
            list_group = group[0] + group[1]
        elif type(group) == dict :
            list_group = [i for i in group.keys()]
        else :
            raise ValueError(f"group type : {type(group)} not an dict or a tuple.")

        print("list_group : ", list_group)
        set_sample = set(self.list_samples)
        set_group = set(list_group)
        missing_elements = set_sample - set_group

        if missing_elements:
            raise ValueError(f"The following elements from set_sample are not present in set_group: {missing_elements}")

    def binary_table(self, snarls, binary_groups) :
        self.check_pheno_group(binary_groups)
        res_table = []
        for snarl in snarls :
            df = self.create_binary_table(binary_groups, snarl)
            res_table.append(df)

        return res_table
    
    def quantitative_table(self, snarls, pheno) :
        self.check_pheno_group(pheno)
        res_table = []
        for snarl in snarls :
            df = self.create_quantitative_table(snarl)
            p_value = self.linear_regression(df, pheno)
            res_table.append((df, p_value))

        return res_table

    def create_quantitative_table(self, snarl) -> pd.DataFrame :
        column_headers = snarl[1]
        row_headers = self.matrix.get_row_header()
        column_headers_header = self.list_samples
        genotypes = []

        # Iterate over each path_snarl in column_headers
        for idx_g, path_snarl in enumerate(column_headers):
            idx_srr_save = list(range(len(column_headers_header)))  # Initialize with all indices
            decomposed_snarl = self.decompose_string(path_snarl)

            # print("path_snarl : ", path_snarl)
            # print("idx_g : ", idx_g)
            # print("idx_srr_save : ", idx_srr_save)
            # print("decomposed_snarl : ", decomposed_snarl)
            genotype = []

            # Iterate over each snarl in the decomposed_snarl
            for snarl in decomposed_snarl:
                
                # print("snarl : ", snarl)
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
                    
            genotype = [1 if idx_header in idx_srr_save else 0 for idx_header in range(len(column_headers_header))]
            genotypes.append(genotype)

        # Transposing the matrix
        transposed_genotypes = [list(row) for row in zip(*genotypes)]
        df = pd.DataFrame(transposed_genotypes, index=column_headers_header, columns=column_headers)
        print(df)
        return df
 
    def create_binary_table(self, groups, snarl) -> pd.DataFrame :
        column_headers = snarl[1]
        row_headers = self.matrix.get_row_header()
        column_headers_header = self.matrix.get_list_samples()
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
        print(df)
        return df

    def linear_regression(self, df, pheno) :

        df = df.astype(int)

        print("pheno : ", pheno)
        print('df : ', df)

        df['Target'] = df.index.map(pheno)
        print('df 2 : ', df)

        # Separate features (X) and target (y)
        X = df.drop('Target', axis=1)
        y = df['Target']
        print("X : ", X)
        print("Y : ", y)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

        model = LogisticRegression()

        # Train the model
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)

        # Evaluate the model's performance
        mse = mean_squared_error(y_test, y_pred)
        r2 = r2_score(y_test, y_pred)
        print(f"Mean Squared Error: {mse}")
        print(f"R^2 Score: {r2}")
        return r2

def parse_group_file(groupe_file) :
    group_0 = []
    group_1 = []
    
    with open(groupe_file, 'r') as file:
        for line in file:
            if line.strip():
                g0, g1 = line.split()
                group_0.append(g0)
                group_1.append(g1)
    
    return group_0, group_1

def parse_pheno_file(phenotype_file) -> dict :
    parsed_data = {}

    with open(phenotype_file, 'r') as file:
        for line in file:
            # Strip any leading/trailing whitespace and split the line by whitespace
            parts = line.strip().split()
            # Add the resulting list to the parsed_data list
            parsed_data[parts[0]] = float(parts[1])

    return parsed_data

def parse_snarl_path_file(path_file) -> list :
    path_list = []
    
    with open(path_file, 'r') as file:
        for line in file:
            if line.strip():  # Skip empty lines
                snarl, aT = line.split()
                path_list.append([snarl, aT.split(',')])
    
    return path_list

def check_format_vcf_file(value):
    """
    Custom function to check if the provided file path is a valid VCF file.
    """
    if not os.path.isfile(value):
        raise argparse.ArgumentTypeError(f"The file {value} does not exist.")

    if not value.lower().endswith('.vcf') or value.lower().endswith('.vcf.gz'):
        raise argparse.ArgumentTypeError(f"The file {value} is not a valid VCF file. It must have a .vcf extension or .vcf.gz.")
    
    return value

def check_format_group_snarl(value):
    """
    Custom function to check if the provided file path is a valid group/snarl file.
    """
    if not os.path.isfile(value):
        raise argparse.ArgumentTypeError(f"The file {value} does not exist.")
    if not value.lower().endswith('.txt') or value.lower().endswith('.tsv'):
        raise argparse.ArgumentTypeError(f"The file {value} is not a valid group/snarl file. It must have a .txt extension or .tsv.")
    
    return value

if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description="Parse and analyse snarl from vcf file")
    parser.add_argument("vcf_path", type=check_format_vcf_file, help="Path to the vcf file")
    parser.add_argument("snarl", type=check_format_group_snarl, help="Path to the snarl file that containt snarl and aT")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--binary", type=check_format_group_snarl, help="Path to the binary group file")
    group.add_argument("-q", "--quantitative", type=check_format_group_snarl, help="Path to the quantitative phenotype file")
    args = parser.parse_args()

    start = time.time()
    vcf_object = Snarl(args.vcf_path)
    vcf_object.initialise_matrix()
    snarl = parse_snarl_path_file(args.snarl)

    if args.binary:
        binary_group = parse_group_file(args.binary)
        vcf_object.binary_table(snarl, binary_group)

    if args.quantitative:
        quantitative = parse_pheno_file(args.quantitative)
        vcf_object.quantitative_table(snarl, quantitative)

    print(f"Time : {time.time() - start} s")

