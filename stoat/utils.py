from cyvcf2 import VCF # type: ignore
import argparse
import numpy as np # type: ignore
import pandas as pd # type: ignore
from collections import defaultdict
import os

# ----------------- Parse functions -----------------

def parsing_samples_vcf(vcf_path:str) -> list :
    try :
        return VCF(vcf_path).samples
    except :
        raise f"Error : VCF parsing error, verify {vcf_path} format"

def parse_covariate_file(filepath:str) -> dict:

    covariates = pd.read_csv(filepath)
    return covariates.set_index("ID").to_dict(orient="index")

def parse_pheno_binary_file(group_file:str) -> dict:

    df = pd.read_csv(group_file, sep='\t')
    
    # Get unique phenotypes
    unique_phenotypes = df['PHENO'].unique()
    
    # Check if there are exactly two phenotypes
    if len(unique_phenotypes) != 2:
        raise ValueError(f"Expected exactly 2 unique phenotypes, but found {len(unique_phenotypes)}: {unique_phenotypes}")
    
    # Map phenotypes to binary values (0 and 1)
    phenotype_mapping = {phenotype: idx for idx, phenotype in enumerate(sorted(unique_phenotypes))}
    
    # Create a dictionary with sample names (IID) as keys and binary phenotypes as values
    binary_pheno = df.set_index('IID')['PHENO'].map(phenotype_mapping).to_dict()
    return binary_pheno

def parse_pheno_quantitatif_file(file_path:str) -> dict:

    df = pd.read_csv(file_path, sep='\t')

    # Extract the IID (second column) and PHENO (third column) and convert PHENO to float
    quantitative_pheno = dict(zip(df['IID'], df['PHENO']))
    return quantitative_pheno

def parse_snarl_path_file(path_file:str) -> tuple[dict, int]:
    
    # Initialize an empty dictionary for the snarl paths
    snarl_paths = defaultdict(list)
    snarl_number_analysis = 0

    # Read the file into a pandas DataFrame
    df = pd.read_csv(path_file, sep='\t', dtype=str)
    df['paths'] = df['paths'].str.split(',')

    # Create the dictionary with snarl as key and paths as values
    for snarl, paths in zip(df['snarl'], df['paths']):
        snarl_paths[snarl].extend(paths)
        snarl_number_analysis += 1

    return snarl_paths, snarl_number_analysis

def parse_covariate_file(covar_path:str) -> dict:

    covariates = pd.read_csv(covar_path)

    # Convert to dictionary with ID as the key and the rest of the row as the value
    return covariates.set_index("ID").to_dict(orient="index")

def parse_plink_grm(prefix: str) -> pd.DataFrame:
    """ Parse PLINK GRM binary files and return the kinship matrix as a pandas DataFrame. """

    # File paths
    grm_bin_file = f"{prefix}.grm.bin"
    grm_id_file = f"{prefix}.grm.id"
    
    # Read IDs (FIDs and IIDs)
    with open(grm_id_file, 'r') as f:
        ids = [line.strip().split() for line in f.readlines()]
        ids = [f"{fid}_{iid}" for fid, iid in ids]  # Combine FID and IID as unique individual IDs
    
    n = len(ids)  # Number of individuals
    
    # Read GRM values
    with open(grm_bin_file, 'rb') as f:
        grm_values = np.frombuffer(f.read(), dtype=np.float64)
    
    # Reconstruct the full GRM matrix (symmetric matrix)
    grm_matrix = np.zeros((n, n))
    idx = 0
    for i in range(n):
        for j in range(i + 1):
            grm_matrix[i, j] = grm_values[idx]
            grm_matrix[j, i] = grm_values[idx]  # Ensure matrix is symmetric
            idx += 1
    
    # Create a pandas DataFrame for easy manipulation
    kinship_matrix = pd.DataFrame(grm_matrix, index=ids, columns=ids)
    return kinship_matrix, ids

# ----------------- Check functions -----------------

def check_file(file_path: str) -> str:
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"Error: File '{file_path}' not found or is not a valid file.")
    return file_path

def check_kinship_prefix(prefix_path: str) -> str:
    """ Function to check if the provided file path is a valid PLINK kinship file. """

    # Check if the other files exist
    check_file(f"{prefix_path}.grm.bin")
    check_file(f"{prefix_path}.grm.id")
    
    return prefix_path

def check_matching(elements: dict, list_samples: list, file_name: str) -> None:
    """Check if all sample names in the pheno file are matching with VCF sample names, else return error."""
    
    set_samples = set(list_samples)
    set_elements = set(elements.keys())
    missing_elements = set_samples - set_elements

    if missing_elements:
        raise ValueError(f"The following sample names from VCF are not present in {file_name} file: {missing_elements}")

def check_format_list_path(file_path:str) -> str:
    """
    Function to check if the provided file path is a valid list path file.
    """
    check_file(file_path)

    with open(file_path, 'r') as file:
        # Read and validate the header
        first_line = file.readline().strip()
        expected_header = 'snarl\tpaths\ttype'
        if first_line != expected_header:
            raise argparse.ArgumentTypeError(
                f"The file must start with the following header: '{expected_header}' and be split by tabulation"
            )
        
        # Validate all other lines
        for line_number, line in enumerate(file, start=2):  # Start at 2 for line number after header
            columns = line.strip().split('\t')
            if len(columns) != 3:
                raise argparse.ArgumentTypeError(
                    f"Line {line_number} must contain exactly 3 columns, but {len(columns)} columns were found."
                )
            if not all(isinstance(col, str) and col.strip() for col in columns):
                raise argparse.ArgumentTypeError(
                    f"Line {line_number} contains empty or non-string values: {columns}"
                )

    return file_path

def check_format_vcf_file(file_path:str) -> str:
    """
    Function to check if the provided file path is a valid VCF file.
    """
    check_file(file_path)

    if not file_path.lower().endswith('.vcf') and not file_path.lower().endswith('.vcf.gz'):
        raise argparse.ArgumentTypeError(f"The file {file_path} is not a valid VCF file. It must have a .vcf extension or .vcf.gz.")
    return file_path

def check_format_pheno(file_path:str) -> str:
    
    check_file(file_path)
    
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
        header = first_line.split('\t')
        expected_header = ['FID', 'IID', 'PHENO']
        if header != expected_header:
            raise argparse.ArgumentTypeError(
                f"The file must contain the following headers: {expected_header} and be split by tabulation"
            )
        
        # Validate all other lines
        for line_number, line in enumerate(file, start=2):  # Start at 2 for line number after header
            columns = line.strip().split('\t')
            if len(columns) != 3:
                raise argparse.ArgumentTypeError(
                    f"Line {line_number} must contain exactly 3 columns, but {len(columns)} columns were found."
                )
            try:
                float(columns[2])  # Check if the third column is a float or can be converted to one
            except ValueError:
                raise argparse.ArgumentTypeError(
                    f"Line {line_number} contains a non-numeric value in the last column: {columns[2]}"
                )

    return file_path

def check_covariate_file(file_path:str) -> str:
    
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
