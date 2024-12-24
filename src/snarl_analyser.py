import argparse
import utils
from cyvcf2 import VCF
from typing import List
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
# from limix.stats import logisticMixedModel
# from limix.stats import scan
import time
import os
 
class Matrix :
    def __init__(self, default_row_number=1_000_000, column_number=2):
        self.default_row_number = default_row_number 
        self.matrix = np.zeros((default_row_number, column_number), dtype=bool)
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

class SnarlProcessor:
    def __init__(self, vcf_path: str, list_samples: list):
        self.list_samples = list_samples
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

    def binary_table(self, snarls, binary_groups, covar=None, gaf=False, output="output/binary_output.tsv"):
        """
        Generate a binary table with statistical results and write to a file.
        """
        headers = ('CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tP_FISHER\tP_CHI2\tTOTAL_SUM\tMIN_ROW_INDEX\tNUM_COLUM\tINTER_GROUP\tAVERAGE\tGROUP_PATHS\n')
        with open(output, 'wb') as outf:
            outf.write(headers.encode('utf-8'))
                
            for snarl, list_snarl in snarls.items():
                # Create the binary table, considering covariates if provided
                df = self.create_binary_table(binary_groups, list_snarl)

                # Perform statistical tests and compute descriptive statistics
                fisher_p_value, chi2_p_value, total_sum, min_sample, numb_colum, inter_group, average, group_paths = self.binary_stat_test(df, gaf)
                chrom = pos = type_var = ref = alt = "NA"
                data = (f'{chrom}\t{pos}\t{snarl}\t{type_var}\t{ref}\t{alt}\t{fisher_p_value}\t'
                        f'{chi2_p_value}\t{total_sum}\t{min_sample}\t{numb_colum}\t{inter_group}\t{average}\t{group_paths}\n')
                outf.write(data.encode('utf-8'))

    def quantitative_table(self, snarls, quantitative, covar=None, output="output/quantitative_output.tsv") :

        with open(output, 'wb') as outf:
            headers = 'CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tBETA\tSE\tP\n'
            outf.write(headers.encode('utf-8'))
            for snarl, list_snarl in snarls.items() :
                df = self.create_quantitative_table(list_snarl)
                beta, se, pvalue = self.linear_regression(df, quantitative)
                chrom = pos = type_var = ref = alt = "NA"
                data = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, snarl, type_var, ref, alt, beta, se, pvalue)
                outf.write(data.encode('utf-8'))

    def identify_correct_path(self, decomposed_snarl: list, idx_srr_save: list) -> list:
        """
        Return a list of column indices where all specific elements of this column in the matrix are 1.
        """

        rows_to_check = np.array([], dtype=int)

        # Print the decomposed_snarl and row_headers_dict
        for snarl in decomposed_snarl:
            if "*" in snarl:
                continue
            if snarl in self.matrix.get_row_header():
                row_index = self.matrix.get_row_header()[snarl]
                rows_to_check = np.append(rows_to_check, row_index)
            else:
                return []

        # Extract the rows from the matrix using rows_to_check
        extracted_rows = self.matrix.get_matrix()[rows_to_check, :]

        # Check if all elements in the columns are 1 for the specified rows
        columns_all_ones = np.all(extracted_rows == 1, axis=0)

        # Find the column indices where all elements are 1
        idx_srr_save = np.where(columns_all_ones)[0].tolist()

        return idx_srr_save

    def create_binary_table(self, groups, list_path_snarl) -> pd.DataFrame:
        """Generates a binary table DataFrame indicating the presence of snarl paths in given groups based on matrix data"""
        
        length_column_headers = len(list_path_snarl)

        # Initialize g0 and g1 with zeros, corresponding to the length of column_headers
        g0 = [0] * length_column_headers
        g1 = [0] * length_column_headers

        # Iterate over each path_snarl in column_headers
        for idx_g, path_snarl in enumerate(list_path_snarl):
            idx_srr_save = list(range(len(self.list_samples)))
            decomposed_snarl = self.decompose_string(path_snarl)
            idx_srr_save = self.identify_correct_path(decomposed_snarl, idx_srr_save)

            # Count occurrences in g0 and g1 based on the updated idx_srr_save
            for idx in idx_srr_save:
                srr = self.list_samples[idx // 2]

                if srr in groups[0]:
                    g0[idx_g] += 1

                if srr in groups[1]:
                    g1[idx_g] += 1

        # Create and return the DataFrame
        df = pd.DataFrame([g0, g1], index=['G0', 'G1'], columns=list_path_snarl)
        return df

    def create_quantitative_table(self, column_headers: list) -> pd.DataFrame:
        length_sample = len(self.list_samples)

        # Initialize a zero matrix for genotypes with shape (length_sample, len(column_headers))
        genotypes = np.zeros((length_sample, len(column_headers)), dtype=int)

        # Iterate over each path_snarl and fill in the matrix
        for col_idx, path_snarl in enumerate(column_headers):
            decomposed_snarl = self.decompose_string(path_snarl)
            idx_srr_save = self.identify_correct_path(decomposed_snarl, list(range(length_sample)))

            for idx in idx_srr_save:
                srr_idx = idx // 2  # Convert index to the appropriate sample index
                genotypes[srr_idx, col_idx] += 1

        df = pd.DataFrame(genotypes, index=self.list_samples, columns=column_headers)
        return df
    
    def linear_regression(self, df, pheno : dict, covar=None) -> tuple :

        df = df.astype(int)
        df['Target'] = df.index.map(pheno)
        x = df.drop('Target', axis=1)
        y = df['Target']

        x_with_const = sm.add_constant(x)
        result = sm.OLS(y, x_with_const).fit()
        beta = result.params  # Beta coefficients
        standard_errors = result.bse  # Standard errors
        formatted_p_value = f"{result.f_pvalue:.4e}"

        return beta, standard_errors, formatted_p_value

    # # Linear Mixed Model
    # def LMM_quantitatif(self, kinship_matrix, covar: dict, pheno: dict) -> tuple:
    #     """
    #     Perform Linear Mixed Model (LMM) for quantitative phenotype data.
    #     """

    #     # Ensure the covariate matrix is in a DataFrame and map the covariates and phenotype correctly
    #     covar_df = pd.DataFrame(covar)
    #     covar_df['Target'] = covar_df.index.map(pheno)  # Map phenotype values
        
    #     # Extract dependent (y) and independent variables (x)
    #     y = covar_df['Target']
    #     x = covar_df.drop('Target', axis=1)  # Remove target from covariates
    #     x = sm.add_constant(x)               # Add constant for intercept term
        
    #     # Perform Linear Mixed Model scan (assuming a function like scan)
    #     results = scan(y=y, K=kinship_matrix, covariates=x)
        
    #     # Extract metrics from the results object (p-value, beta, beta_se, log-likelihood, heritability)
    #     p_value = round(results.stats["pv"], 4)  # P-values for each covariate
    #     beta = results.stats["beta"]            # Effect sizes (coefficients for covariates)
    #     beta_se = results.stats["beta_se"]      # Standard errors for effect sizes
    #     ll = results.stats["ll"]                # Log-likelihood of the model
    #     heritability = results.stats["h2"]      # Heritability estimate (proportion of variance explained by GRM)
    #     return p_value, beta, beta_se, ll, heritability

    def chi2_test(self, df) -> float:
        """Calculate p_value from list of dataframe using chi-2 test"""

        # Check if dataframe has at least 2 columns and more than 0 counts in every cell
        if df.shape[1] >= 2 and np.all(df.sum(axis=0)) and np.all(df.sum(axis=1)):
            try:
                # Perform Chi-Square test
                p_value = chi2_contingency(df)[1] # from scipy.stats import chi2_contingency
                p_value = f"{p_value:.4e}"

            except ValueError as e:
                p_value = "Error"
        else:
            p_value = "N/A"

        return p_value

    def fisher_test(self, df) -> float : 
        """Calcul p_value using fisher exact test"""

        try:
            p_value = fisher_exact(df)[1] # from scipy.stats import fisher_exact
            p_value = f"{p_value:.4e}"

        except ValueError as e: 
            p_value = 'N/A'
        
        return p_value

    # # Logistic Mixed Model
    # def LMM_binary(df, pheno, covar):
    #     """
    #     Perform Logistic Mixed Model on a binary phenotype.
    #     """
    #     # Ensure phenotype mapping
    #     if not all(ind in df.index for ind in pheno.keys()):
    #         raise ValueError("Some individuals in the phenotype file do not match the genotype file.")

    #     # Map phenotype to df
    #     df['Target'] = df.index.map(pheno)
    #     y = df['Target'].values
    #     X = covar[df.index].values  # Covariates should match the index of genotype data
    #     K = np.corrcoef(df.T)       # This is a placeholder for the kinship matrix (should be computed properly in practice)
    #     lmm = logisticMixedModel(y=y, K=K)
        
    #     # Fit the model with covariates
    #     lmm.fit(X)
    #     beta = lmm.beta                 # Effect sizes (log-odds for covariates)
    #     p_value = round(lmm.pv, 4)      # P-values for fixed effects
    #     vcomp = lmm.vcomp               # Variance components (relatedness)

    #     return p_value, beta, vcomp
    
    def format_group_paths(self, df) : 
        result = []
        for column in df.columns:
            column_values = [f"{df.loc[group, column]}" for group in df.index]
            result.append(":".join(column_values))

        # Join all results with ", "
        final_str = ",".join(result)
        return final_str
    
    def binary_stat_test(self, df, gaf) :
        """ Perform statistical tests and calculate descriptive statistics on a binary data frame. """
        fisher_p_value = self.fisher_test(df)
        chi2_p_value = self.chi2_test(df)

        total_sum = int(df.values.sum()) # Calculate the total sum of all values in the DataFrame
        inter_group = int(df.min().sum()) # Calculate the sum of the minimum values from each column
        numb_colum = df.shape[1] # Get the total number of columns in the DataFrame
        average = float(total_sum / numb_colum) # Calculate the average of the total sum divided by the number of columns
        row_sums = df.sum(axis=1) # Calculate the sum of values for each row in the DataFrame
        min_row_index = row_sums.min() # Find the minimum sum among the row sums
        group_paths = self.format_group_paths(df) if gaf else ""

        return fisher_p_value, chi2_p_value, total_sum, min_row_index, numb_colum, inter_group, average, group_paths
    
if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description="Parse and analyse snarl from vcf file")
    parser.add_argument("vcf_path", type=utils.check_format_vcf_file, help="Path to the vcf file (.vcf or .vcf.gz)")
    parser.add_argument("snarl", type=utils.check_format_list_path, help="Path to the snarl file that containt snarl and aT (.txt or .tsv)")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--binary", type=utils.check_format_pheno, help="Path to the binary phenotype file (.txt or .tsv)")
    group.add_argument("-q", "--quantitative", type=utils.check_format_pheno, help="Path to the quantitative phenotype file (.txt or .tsv)")
    parser.add_argument("-c", "--covariate", type=utils.check_covariate_file, required=False, help="Path to the covariate file (.txt or .tsv)")
    parser.add_argument("-o", "--output", type=str, required=False, help="Path to the output file")
    args = parser.parse_args()

    start = time.time()
    list_samples = utils.parsing_samples_vcf(args.vcf_path)
    vcf_object = SnarlProcessor(args.vcf_path, list_samples)
    vcf_object.fill_matrix()
    print(f"Time Matrix : {time.time() - start} s")

    start = time.time()
    snarl = utils.parse_snarl_path_file(args.snarl)[0]
    covar = utils.parse_covariate_file(args.covariate) if args.covariate else ""

    output_dir = args.output or "output"    
    os.makedirs(output_dir, exist_ok=True)
    output = os.path.join(output_dir, "stoat.assoc.tsv")

    if args.binary:
        binary_group = utils.parse_pheno_binary_file(args.binary)
        vcf_object.binary_table(snarl, binary_group, covar, output)

    # python3 src/snarl_analyser.py tests/other_files/small_vcf.vcf tests/other_files/list_snarl_short.txt -b tests/other_files/group.txt
    # python3 src/snarl_analyser.py ../snarl_data/fly.merged.vcf output/test_list_snarl.tsv -b ../snarl_data/group.txt
    # python3 src/snarl_analyser.py ../snarl_data/simulation_1000vars_100samps/calls/merged_output.vcf ../snarl_data/simulation_1000vars_100samps/pg.snarl_netgraph.paths.tsv -b ../snarl_data/simulation_1000vars_100samps/group.txt -o output/simulation_binary.tsv

    if args.quantitative:
        quantitative = utils.parse_pheno_quantitatif_file(args.quantitative)
        vcf_object.quantitative_table(snarl, quantitative, covar, output)

    # python3 src/snarl_analyser.py tests/other_files/small_vcf.vcf tests/other_files/list_snarl_short.txt -q tests/other_files/pheno.txt

    print(f"Time P-value: {time.time() - start} s")


