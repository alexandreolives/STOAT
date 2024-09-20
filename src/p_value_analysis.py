import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 

def plot_p_value_distribution_binary(file_path, output_file_dist):
    df = pd.read_csv(file_path, sep='\t')

    p_values_f = df.iloc[:, 0]
    p_values_c = df.iloc[:, 1]

    plt.figure(figsize=(10, 6))
    plt.hist(p_values_f, bins=50, color='skyblue', edgecolor='black', alpha=0.7, label='P-value (Fisher)')
    plt.hist(p_values_c, bins=50, color='orange', edgecolor='black', alpha=0.7, label='P-value (Chi2)')
    plt.title('P-value Distribution', fontsize=16)
    plt.xlabel('P-value', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.legend(loc='upper left')
    plt.savefig(output_file_dist, format='png', dpi=300)

def plot_p_value_distribution_quantitative(file_path, output_file_dist):

    df = pd.read_csv(file_path, sep='\t')
    p_values = df.iloc[:, 1]

    plt.figure(figsize=(10, 6))
    plt.hist(p_values, bins=50, color='skyblue', edgecolor='black', alpha=0.7, label='P-value')
    plt.title('P-value Distribution', fontsize=16)
    plt.xlabel('P-value', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.legend(loc='upper left')
    plt.savefig(output_file_dist, format='png', dpi=300)

def plot_manhattan_binary(file_path, output_file_manh):
    
    df = pd.read_csv(file_path, sep='\t')
    
    # Convert 'N/A' to NaN and then to numeric, while filtering out NaN values
    df['Snarl'] = pd.to_numeric(df['Snarl'], errors='coerce')
    df['P_value_Fisher'] = pd.to_numeric(df['P_value_Fisher'], errors='coerce')
    df = df.dropna(subset=['Snarl', 'P_value_Fisher'])

    # Convert p-values to -log10(p-values) for plotting
    df['-log10(P_value_Fisher)'] = -np.log10(df['Snarl'])
    df['-log10(P_value_Chi2)'] = -np.log10(df['P_value_Fisher'])
    df['index'] = np.arange(len(df))
    
    plt.figure(figsize=(12, 6))
    plt.plot(df['index'], df['-log10(P_value_Fisher)'], 'o', color='skyblue', alpha=0.7, label='P_value_Fisher')
    plt.plot(df['index'], df['-log10(P_value_Chi2)'], 'o', color='orange', alpha=0.7, label='P_value_Chi2')
    
    significance_threshold_1 = -np.log10(0.00001)
    significance_threshold_2 = -np.log10(0.00000001)

    plt.axhline(y=significance_threshold_1, color='red', linestyle='dashed', label='p=0.00001')
    plt.axhline(y=significance_threshold_2, color='green', linestyle='dashed', label='p=0.00000001')

    plt.title('Manhattan Plot', fontsize=16)
    plt.xlabel('Index', fontsize=14)
    plt.ylabel('-log10(P-value)', fontsize=14)
    plt.legend(loc='upper right')
    plt.savefig(output_file_manh, format='png', dpi=300)

def plot_manhattan_quantitative(file_path, output_file_manh):
    
    df = pd.read_csv(file_path, sep='\t')
    
    df['P_value'] = pd.to_numeric(df['P_value'], errors='coerce')
    df = df.dropna(subset=['P_value'])
    df['-log10(P_value)'] = -np.log10(df['P_value'])
    df['index'] = np.arange(len(df))

    plt.figure(figsize=(12, 6))
    plt.plot(df['index'], df['-log10(P_value)'], 'o', color='skyblue', alpha=0.7, label='P-value')
    significance_threshold = -np.log10(0.00001)

    plt.axhline(y=significance_threshold, color='red', linestyle='--', label='p=0.00001')
    plt.title('Manhattan Plot of P-values', fontsize=16)
    plt.xlabel('Index', fontsize=14)
    plt.ylabel('-log10(P-value)', fontsize=14)
    plt.legend(loc='upper right')
    plt.savefig(output_file_manh, format='png', dpi=300)

# file_path = 'output/pgtest_100_binary.tsv'
# output_file_dist = "output/pgtest_100_distribution_plot_binary.png"
# output_file_manh = "output/pgtest_100_manhattan_plot_binary.png"
# plot_p_value_distribution_binary(file_path, output_file_dist)
# plot_manhattan_binary(file_path, output_file_manh)

file_path = 'output/simulation_quantitative.tsv'
output_file_dist = "output/pgtest_1000_distribution_plot_quantitative.png"
output_file_manh = "output/pgtest_1000_manhattan_plot_quantitative.png"
plot_p_value_distribution_quantitative(file_path, output_file_dist)
plot_manhattan_quantitative(file_path, output_file_manh)

# file_path = 'output/droso_quantitative.tsv'
# output_file_dist = "output/distribution_plot_quantitative.png"
# output_file_manh = "output/manhattan_plot_quantitative.png"
# plot_p_value_distribution_quantitative(file_path, output_file_dist)
# plot_manhattan_quantitative(file_path, output_file_manh)
