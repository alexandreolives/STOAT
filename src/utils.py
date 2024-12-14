from cyvcf2 import VCF
import argparse
import pandas as pd
from collections import defaultdict
import os

def parsing_samples_vcf(vcf_path) :
    try :
        return VCF(vcf_path).samples
    except :
        raise f"Error : VCF parsing error, verify {vcf_path} format"


def parse_covariate_file(filepath):

    covariates = pd.read_csv(filepath)
    return covariates.set_index("ID").to_dict(orient="index")

def parse_pheno_binary_file(group_file : str):

    df = pd.read_csv(group_file, sep='\t')
    
    # Ensure binary phenotype is valid
    if not set(df['PHENO'].dropna()).issubset({0, 1}):
        raise ValueError("The 'PHENO' column must contain only binary values (0 or 1).")

    # Create dictionaries for group 0 and group 1
    group_0 = {sample: 0 for sample in df[df['PHENO'] == 0]['IID']}
    group_1 = {sample: 1 for sample in df[df['PHENO'] == 1]['IID']}
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
    
def check_format_pheno(file_path : str) -> str :

    check_file(file_path)
    
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
    
    header = first_line.split('\t')
    expected_header = ['FID', 'IID', 'PHENO']
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