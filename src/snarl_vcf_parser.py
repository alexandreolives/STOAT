import argparse
from cyvcf2 import VCF
from typing import List
import numpy as np
import pandas as pd
import statsmodels.api as sm
from collections import defaultdict
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
from limix.stats import logisticMixedModel
import os
import time

def parsing_samples_vcf(vcf_path) :
    try :
        return VCF(vcf_path).samples
    except :
        raise f"Error : VCF parsing error, verify {vcf_path} format"

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

    def __str__(self) :
        return f"{self.matrix[0]} \n" \
               f"{self.matrix[1]} \n" \
               f"{self.matrix[2]} \n" \
               f"{self.matrix[3]} \n" \
               f"{self.matrix[4]} \n"

class SnarlProcessor:
    def __init__(self, vcf_path: str, list_samples):
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

    def binary_table(self, snarls, binary_groups, covar=None, output="output/binary_output.tsv"):
        """
        Generate a binary table with statistical results and write to a file.
        """
        headers = ('CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tP_FISHER\tP_CHI2\tTOTAL_SUM\t'
                'MIN_ROW_INDEX\tINTER_GROUP\tAVERAGE\n')
        
        with open(output, 'wb') as outf:
            outf.write(headers.encode('utf-8'))
            
            for snarl, list_snarl in snarls.items():
                # Create the binary table, considering covariates if provided
                if covar:
                    df = self.create_binary_table(binary_groups, list_snarl, covar)
                else:
                    df = self.create_binary_table(binary_groups, list_snarl)
                
                # Perform statistical tests and compute descriptive statistics
                fisher_p_value, chi2_p_value, total_sum, min_row_index, inter_group, average = self.binary_stat_test(df)
                chrom = pos = type_var = ref = alt = "NA"
                data = (f'{chrom}\t{pos}\t{snarl}\t{type_var}\t{ref}\t{alt}\t{fisher_p_value}\t'
                        f'{chi2_p_value}\t{total_sum}\t{min_row_index}\t{inter_group}\t{average}\n')
                outf.write(data.encode('utf-8'))

    def quantitative_table(self, snarls, quantitative, covar=None, output="output/quantitative_output.tsv") :

        with open(output, 'wb') as outf:
            headers = 'CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tP\n'
            outf.write(headers.encode('utf-8'))
            for snarl, list_snarl in snarls.items() :
                df = self.create_quantitative_table(list_snarl)
                pvalue = self.linear_regression(df, quantitative)
                chrom = pos = type_var = ref = alt = "NA"
                data = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, snarl, type_var, ref, alt, pvalue)
                outf.write(data.encode('utf-8'))

    def identify_correct_path(self, decomposed_snarl: list, row_headers_dict: dict, idx_srr_save: list) -> list:
        """
        Return a list of column indices where all specific elements of this column in the matrix are 1.
        """
        matrix = self.matrix.get_matrix()  # Assuming this returns the matrix as a numpy array
        rows_to_check = np.array([], dtype=int)

        # Print the decomposed_snarl and row_headers_dict

        for snarl in decomposed_snarl:
            if "*" in snarl:
                continue
            if snarl in row_headers_dict:
                row_index = row_headers_dict[snarl]
                rows_to_check = np.append(rows_to_check, row_index)
            else:
                return []

        # Extract the rows from the matrix using rows_to_check
        extracted_rows = matrix[rows_to_check, :]

        # Check if all elements in the columns are 1 for the specified rows
        columns_all_ones = np.all(extracted_rows == 1, axis=0)

        # Find the column indices where all elements are 1
        idx_srr_save = np.where(columns_all_ones)[0].tolist()

        return idx_srr_save

    def create_binary_table(self, groups, list_path_snarl, covar=None) -> pd.DataFrame:
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
            for idx in idx_srr_save:
                srr = list_samples[idx // 2]

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

        # Fit the regression model
        x_with_const = sm.add_constant(x)
        result = sm.OLS(y, x_with_const).fit()
        return result.f_pvalue
    
    def linear_regression(self, df, pheno : dict, covar=None) -> float :
        
        df = df.astype(int)
        df['Target'] = df.index.map(pheno)
        x = df.drop('Target', axis=1)
        y = df['Target']
        pval = round(self.sm_ols(x, y), 4)
        return pval

    def compute_kinship_matrix(self, df):
        """
        Compute the kinship matrix for a given genotypic dataset.
        """
        # Compute allele frequencies for each SNP
        allele_frequencies = df.mean() / 2  # Average genotype value, assuming 0, 1, 2 coding.

        # Center the genotypes by subtracting allele frequencies
        centered_genotypes = df - allele_frequencies

        # The kinship coefficient between individuals i and j is the sum of the product of centered genotypes across SNPs, divided by 2
        kinship_matrix = centered_genotypes.T.dot(centered_genotypes) / (2 * df.shape[1])
        return kinship_matrix

    # Linear Mixed Model
    def LMM_quantitatif(self, kinship_matrix, covar: dict, pheno: dict) -> tuple:
        """
        Perform Linear Mixed Model (LMM) for quantitative phenotype data.
        """

        # Ensure the covariate matrix is in a DataFrame and map the covariates and phenotype correctly
        covar_df = pd.DataFrame(covar)
        covar_df['Target'] = covar_df.index.map(pheno)  # Map phenotype values
        
        # Extract dependent (y) and independent variables (x)
        y = covar_df['Target']
        x = covar_df.drop('Target', axis=1)  # Remove target from covariates
        x = sm.add_constant(x)  # Add constant for intercept term
        
        # Perform Linear Mixed Model scan (assuming a function like scan)
        results = self.scan(y=y, K=kinship_matrix, covariates=x)  # Assuming this is a method that fits the model
        
        # Extract metrics from the results object (p-value, beta, beta_se, log-likelihood, heritability)
        p_value = round(results.stats["pv"],4)  # P-values for each covariate
        beta = results.stats["beta"]  # Effect sizes (coefficients for covariates)
        beta_se = results.stats["beta_se"]  # Standard errors for effect sizes
        ll = results.stats["ll"]  # Log-likelihood of the model
        heritability = results.stats["h2"]  # Heritability estimate (proportion of variance explained by GRM)

        return p_value, beta, beta_se, ll, heritability

    def chi2_test(self, df) -> float:
        """Calculate p_value from list of dataframe using chi-2 test"""

        # Check if dataframe has at least 2 columns and more than 0 counts in every cell
        if df.shape[1] >= 2 and np.all(df.sum(axis=0)) and np.all(df.sum(axis=1)):
            try:
                # Perform Chi-Square test
                p_value = round(chi2_contingency(df)[1], 4) # from scipy.stats import chi2_contingency
            except ValueError as e:
                p_value = "Error"
        else:
            p_value = "N/A"

        return p_value

    def fisher_test(self, df) -> float : 
        """Calcul p_value using fisher exact test"""

        try:
            p_value = round(fisher_exact(df)[1], 4) # from scipy.stats import fisher_exact

        except ValueError as e: 
            p_value = 'N/A'
        
        return p_value

    # Logistic Mixed Model
    def LMM_binary(df, pheno, covar):
        """
        Perform Logistic Mixed Model on a binary phenotype.
        """
        # Ensure phenotype mapping
        if not all(ind in df.index for ind in pheno.keys()):
            raise ValueError("Some individuals in the phenotype file do not match the genotype file.")

        # Map phenotype to df
        df['Target'] = df.index.map(pheno)
        y = df['Target'].values
        X = covar[df.index].values  # Covariates should match the index of genotype data
        K = np.corrcoef(df.T)  # This is a placeholder for the kinship matrix (should be computed properly in practice)
        lmm = logisticMixedModel(y=y, K=K)
        
        # Fit the model with covariates
        lmm.fit(X)
        beta = lmm.beta       # Effect sizes (log-odds for covariates)
        p_value = round(lmm.pv, 4)      # P-values for fixed effects
        vcomp = lmm.vcomp     # Variance components (relatedness)

        return p_value, beta, vcomp
    
    def binary_stat_test(self, df) :
        """ Perform statistical tests and calculate descriptive statistics on a binary data frame. """
        fisher_p_value = self.fisher_test(df)
        chi2_p_value = self.chi2_test(df)
        
        total_sum = int(df.values.sum()) # Calculate the total sum of all elements in the DataFrame
        inter_group = int(df.min().sum()) # Calculate the inter-group value: the sum of the minimum values in each column
        numb_colum = df.shape[1] # Number of columns in the DataFrame
        average = float(total_sum / numb_colum) # Calculate the average value across all elements divided by the number of columns
        row_sums = df.sum(axis=1) # Calculate the sum of each row
        min_row_index = row_sums.min() # Find the minimum row sum across all rows
    
        return fisher_p_value, chi2_p_value, total_sum, min_row_index, inter_group, average

def parse_covariate_file(filepath):
    
    covariates = pd.read_csv(filepath)
    return covariates.set_index("ID").to_dict(orient="index")

def parse_pheno_binary_file(group_file : str) -> tuple:

    df = pd.read_csv(group_file, sep='\t')
    
    # Ensure binary phenotype is valid
    if not set(df['GROUP'].dropna()).issubset({0, 1}):
        raise ValueError("The 'GROUP' column must contain only binary values (0 or 1).")

    # Create dictionaries for group 0 and group 1
    group_0 = {sample: 0 for sample in df[df['GROUP'] == 0]['SAMPLE']}
    group_1 = {sample: 1 for sample in df[df['GROUP'] == 1]['SAMPLE']}
    return group_0, group_1
 
def parse_pheno_quantitatif_file(file_path : str) -> dict:

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

def parse_covariate_file(covar_path : str) -> dict :

    covariates = pd.read_csv(covar_path)

    # Convert to dictionary with ID as the key and the rest of the row as the value
    return covariates.set_index("ID").to_dict(orient="index")

def check_file(file_path) :
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"Error: File '{file_path}' not found or is not a valid file.")
    return file_path

def check_mathing(elements, list_samples, file_name) :
    """Check if all sample name in the pheno file are matching with vcf sample name else return error"""
    
    set_sample = set(list_samples)
    set_elements = set(elements)
    missing_elements = set_sample - set_elements

    if missing_elements :
        raise ValueError(f"The following sample name from vcf are not present in {file_name} file : {missing_elements}")

def check_format_vcf_file(file_path : str) -> str:
    """
    Function to check if the provided file path is a valid VCF file.
    """
    check_file(file_path)

    if not file_path.lower().endswith('.vcf') and not file_path.lower().endswith('.vcf.gz'):
        raise argparse.ArgumentTypeError(f"The file {file_path} is not a valid VCF file. It must have a .vcf extension or .vcf.gz.")
    return file_path

def check_format_group_snarl(file_path : str) -> str :
    """
    Function to check if the provided file path is a valid group file.
    """
    check_file(file_path)
    
    if not file_path.lower().endswith('.txt') and not file_path.lower().endswith('.tsv'):
        raise argparse.ArgumentTypeError(f"The file {file_path} is not a valid group/snarl file. It must have a .txt extension or .tsv.")
    return file_path
    
def check_format_pheno_q(file_path : str) -> str :

    check_file(file_path)
    
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
    
    header = first_line.split('\t')
    expected_header = ['FID', 'IID', 'PHENO']
    if header != expected_header:
        raise argparse.ArgumentTypeError(f"The file must contain the following headers: {expected_header} and be split by tabulation")
    
    return file_path

def check_format_pheno_b(file_path : str) -> str :

    check_file(file_path)

    with open(file_path, 'r') as file:
        first_line = file.readline().strip()

    header = first_line.split('\t')
    expected_header = ['SAMPLE', 'GROUP']
    if header != expected_header:
        raise argparse.ArgumentTypeError(f"The file must contain the following headers: {expected_header} and be split by tabulation")

    return file_path

def check_covariate_file(file_path):
    
    # Check if the file exists
    check_file(file_path)
    
    try:
        covariates = pd.read_csv(file_path, delim_whitespace=True, header=None)
    except Exception as e:
        raise argparse.ArgumentTypeError(f"Error reading the file: {e}")
    
    # Check if the file has at least 5 columns (FID, IID, Age, Sex, PC1)
    if covariates.shape[1] < 5:
        raise argparse.ArgumentTypeError("File must have at least 5 columns: FID, IID, Age, Sex, PC1, PC2")
    
    # Name the columns (PLINK-style file doesn't have headers)
    covariates.columns = ["FID", "IID", "Age", "Sex", "PC1", "PC2"] + [f"PC{i}" for i in range(3, covariates.shape[1])]
    
    # Check for duplicate IDs (FID, IID should be unique)
    if covariates["IID"].duplicated().any():
        raise argparse.ArgumentTypeError("Duplicate IDs found in the file.")
    
    # Check for missing data
    if covariates.isnull().any().any():
        raise argparse.ArgumentTypeError("There are missing values in the file.")
    
    # Check data types
    if not pd.api.types.is_numeric_dtype(covariates["Age"]):
        raise argparse.ArgumentTypeError("Age column should be numeric.")
    
    if covariates["Sex"].dtype != 'object':
        raise argparse.ArgumentTypeError("Sex column should be categorical.")
    
    if not pd.api.types.is_numeric_dtype(covariates["PC1"]) or not pd.api.types.is_numeric_dtype(covariates["PC2"]):
        raise argparse.ArgumentTypeError("PC1 and PC2 columns should be numeric.")
    
    return file_path
    
if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description="Parse and analyse snarl from vcf file")
    parser.add_argument("vcf_path", type=check_format_vcf_file, help="Path to the vcf file (.vcf or .vcf.gz)")
    parser.add_argument("snarl", type=check_format_group_snarl, help="Path to the snarl file that containt snarl and aT (.txt or .tsv)")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--binary", type=check_format_pheno_b, help="Path to the binary group file (.txt or .tsv)")
    group.add_argument("-q", "--quantitative", type=check_format_pheno_q, help="Path to the quantitative phenotype file (.txt or .tsv)")
    parser.add_argument("-c", "--covariate", type=check_covariate_file, required=False, help="Path to the covariate file (.txt or .tsv)")
    parser.add_argument("-o", "--output", type=str, required=False, help="Path to the output file")
    args = parser.parse_args()
    
    start = time.time()
    vcf_object = SnarlProcessor(args.vcf_path)
    vcf_object.fill_matrix()
    print(f"Time Matrix : {time.time() - start} s")

    start = time.time()
    snarl = parse_snarl_path_file(args.snarl)
    covar = parse_covariate_file(args.covariate) if args.covariate else ""

    if args.binary:
        binary_group = parse_pheno_binary_file(args.binary)
        if args.output :
            vcf_object.binary_table(snarl, binary_group, covar, args.output)
        else :
            vcf_object.binary_table(snarl, binary_group, covar)

    # python3 src/snarl_vcf_parser.py test/small_vcf.vcf test/list_snarl_short.txt -b test/group.txt
    # python3 src/snarl_vcf_parser.py ../snarl_data/fly.merged.vcf output/test_list_snarl.tsv -b ../snarl_data/group.txt
    # python3 src/snarl_vcf_parser.py ../snarl_data/simulation_1000vars_100samps/calls/merged_output.vcf ../snarl_data/simulation_1000vars_100samps/pg.snarl_netgraph.paths.tsv -b ../snarl_data/simulation_1000vars_100samps/group.txt -o output/simulation_binary.tsv

    if args.quantitative:
        quantitative = parse_pheno_quantitatif_file(args.quantitative)
        if args.output :
            vcf_object.quantitative_table(snarl, quantitative, covar, args.output)
        else :
            vcf_object.quantitative_table(snarl, quantitative, covar)

    print(f"Time P-value: {time.time() - start} s")


