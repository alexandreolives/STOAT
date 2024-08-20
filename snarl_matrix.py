import io
import argparse
from cyvcf2 import VCF
import numpy as np
import pandas as pd
import os
from typing import List
import time
import pandas as pd
import re 

class Chunk:
    def __init__(self, nb_matrix_row, nb_matrix_column):
        self.data = np.zeros((nb_matrix_row, nb_matrix_column), dtype=bool)

    def add_data(self, idx_snarl, idx_geno):
        self.data[idx_snarl, idx_geno] = 1

    def get_data(self) :
        return self.data
    
    @staticmethod
    def concatenate_matrix(matrix_list, axis=0):
        
        matrix_list = [matrix.get_data() for matrix in matrix_list]

        if not matrix_list:
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
        self.Matrix = Tupple_matrix[0]
        self.row_header = Tupple_matrix[1]
        self.column_header = Tupple_matrix[2]
        
    def __str__(self) :
        return f"           {self.column_header} \n" \
               f"{self.row_header[0]} {self.Matrix[0]}"
        
class Snarl :
    def __init__(self, vcf_path: str, group_file : str, path_file : str) :
        self.paths = self._parse_path_file(path_file)
        self.group = self._parse_group_file(group_file)
        self.matrix = Matrix(self._parse_vcf(vcf_path))

    def print_matrix(self) :
        print(self.matrix)

    def is_valid_vcf(self, file_path: str) -> bool:
        """Check if the file is a valid VCF or VCF.GZ file."""
        return os.path.isfile(file_path) and file_path.endswith(('.vcf', '.vcf.gz'))

    def determine_str(self, s: str, length_s : int, i: int) -> tuple[int, int]:
        """Extract an integer from a string starting at index i."""
        start_idx = i
        while i < length_s and s[i] not in ['>', '<']:
            i += 1
        return i, int(s[start_idx:i])

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

    def push_matrix(self, decomposed_snarl, allele, nb_matrix_row, nb_matrix_column, row_header, Matrix_list, idx_geno, idx_matrix) :
        """Add True to the matrix if snarl are found"""
        for decomposed_snarl_path in decomposed_snarl[allele] :
            if decomposed_snarl_path not in row_header :
                row_header.append(decomposed_snarl_path)
            
            idx_snarl = row_header.index(decomposed_snarl_path)
            
            # Check if current numpy matrix are full 
            if (len(row_header)+1) % nb_matrix_row == 1 :
                # Add new chunck
                Matrix_list.append(Chunk(nb_matrix_row, nb_matrix_column))
                idx_matrix += 1

            Matrix_list[idx_matrix].add_data(idx_snarl%nb_matrix_row, idx_geno) # matrix[idx_snarl, idx_geno] = 1
        
        return Matrix_list, idx_matrix

    def _parse_vcf(self, vcf_path) :
        """Parse vcf file (main function)"""
        if not self.is_valid_vcf(vcf_path) :
            raise ValueError(f"VCF file : {vcf_path} are not on the correct format (.vcf or .vcf.gz)")
        
        vcf_object = VCF(vcf_path)
        list_samples = vcf_object.samples

        # define column and row header for the numpy matrix 
        row_header = []
        column_header = [f"{sample}_{i}" for sample in list_samples for i in range(2)]
        nb_matrix_row = 100000
        nb_matrix_column = len(list_samples) * 2
        Matrix_list = [Chunk(nb_matrix_row, nb_matrix_column)]
        idx_matrix = 0

        # Parse line per line
        for variant in vcf_object :
            genotypes = variant.genotypes # genotypes : [[-1, -1, False], [-1, -1, False], [1, 1, False], [-1, -1, False]]
            snarl_list = variant.INFO.get('AT', '').split(',') # snarl_list :  ['>1272>1273>1274', '>1272>1274']
            decomposed_snarl = self.decompose_snarl(snarl_list) # decomposed_snarl :  [['>1272>1273', '>1273>1274'], ['>1272>1274']]

            idx_geno = 0
            for gt in genotypes:
                allele_1, allele_2 = gt[:2]  # Extract the two allele

                if allele_1 != -1 : # Suppose that allele_2 are always -1
                    Matrix_list, idx_matrix = self.push_matrix(decomposed_snarl, allele_1, nb_matrix_row, nb_matrix_column, row_header, Matrix_list, idx_geno, idx_matrix)
                    Matrix_list, idx_matrix = self.push_matrix(decomposed_snarl, allele_2, nb_matrix_row, nb_matrix_column, row_header, Matrix_list, idx_geno+1, idx_matrix)

                idx_geno += 2 # push the index
            
        return Chunk.concatenate_matrix(Matrix_list), row_header, column_header

    def _parse_group_file(self, groupe_file) :
        group_0 = []
        group_1 = []
        
        with open(groupe_file, 'r') as file:
            for line in file:
                if line.strip():
                    g0, g1 = line.split()
                    group_0.append(g0)
                    group_1.append(g1)
        
        print("group_0, group_1 : ", group_0, group_1)
        return group_0, group_1

    def _parse_path_file(self, path_file) :
        snarl_list = []
        aT_list = []
        
        with open(path_file, 'r') as file:
            for line in file:
                if line.strip():  # Skip empty lines
                    snarl, aT = line.split()
                    snarl_list.append(snarl)
                    aT_list.append(aT.split(','))
        
        print("snarl_list, aT_list : ", snarl_list, aT_list)
        return snarl_list, aT_list

if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description="Parse and analyse snarl from vcf file")
    parser.add_argument("vcf_path", type=str, help="Path to the vcf file")
    parser.add_argument("group", type=str, help="Path to the group file that containt SRR name for each group")
    parser.add_argument("snarl", type=str, help="Path to the snarl file that containt snarl and aT")
    args = parser.parse_args()

    start = time.time()
    vcf_object = Snarl(args.vcf_path, args.group, args.snarl)
    print(vcf_object.print_matrix())
    print(f"Time : {time.time() - start} s")


