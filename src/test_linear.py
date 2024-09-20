import pandas as pd
import statsmodels.api as sm
import argparse
import numpy as np

def sm_ols(x, y):
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    print("x : ", x)
    print("y : ", y)

    result = sm.OLS(y, x).fit()
    print(result.summary())

    list_p_value = []
    for p_value in result.pvalues:
        formatted_p_value = "N/A" if np.isnan(p_value) else p_value
        list_p_value.append(formatted_p_value)

    return list_p_value

def parse_pheno_file(pheno_file_path):
    pheno_data = np.loadtxt(pheno_file_path, delimiter=' ', usecols=[1], skiprows=1)
    pheno_data = pheno_data.reshape(-1, 1)
    return pheno_data

def parse_data_file(data_file_path):
    data = np.loadtxt(data_file_path, dtype=int, delimiter=' ', usecols=[1,2], skiprows=1)
    return data

def main(pheno_file_path, data_file_path):
    y = parse_pheno_file(pheno_file_path)
    x = parse_data_file(data_file_path)

    p_values = sm_ols(x, y)
    print("P-values:", p_values)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse phenotype and data files to perform OLS regression.")
    parser.add_argument('pheno_file', type=str, help='The path to the input phenotype file.')
    parser.add_argument('data_file', type=str, help='The path to the input data file.')
    args = parser.parse_args()
    main(args.pheno_file, args.data_file)
