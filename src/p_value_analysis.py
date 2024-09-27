import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 

def plot_p_value_distribution_binary(file_path, output_file_dist):
    df = pd.read_csv(file_path, sep='\t')

    p_values_f = df['P_value_Fisher']
    p_values_c = df['P_value_Chi2']

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

def plot_manhattan_binary(file_path, output_file_manh, output_snarl="output/significative_snarl.tsv"):
    
    df = pd.read_csv(file_path, sep='\t')
    
    # Convert 'N/A' to NaN and then to numeric, while filtering out NaN values
    df['P_value_Fisher'] = pd.to_numeric(df['P_value_Fisher'], errors='coerce')
    df['P_value_Chi2'] = pd.to_numeric(df['P_value_Chi2'], errors='coerce')
    df = df.dropna(subset=['P_value_Fisher', 'P_value_Chi2'])

    # Convert p-values to -log10(p-values) for plotting
    df['-log10(P_value_Fisher)'] = -np.log10(df['P_value_Fisher'])
    df['-log10(P_value_Chi2)'] = -np.log10(df['P_value_Chi2'])
    df['index'] = np.arange(len(df))
    
    significance_threshold = 0.00001
    tupple_snarl = df[df['P_value'] < significance_threshold][['Snarl', 'P_value_Fisher']].itertuples(index=False, name=None)
    write_significative_snarl(tupple_snarl, output_snarl)

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

def plot_manhattan_quantitative(file_path, output_file_manh, output_snarl="output/significative_snarl.tsv"):
    
    df = pd.read_csv(file_path, sep='\t')
    
    df['P_value'] = pd.to_numeric(df['P_value'], errors='coerce')
    df = df.dropna(subset=['P_value'])
    df['-log10(P_value)'] = -np.log10(df['P_value'])
    df['index'] = np.arange(len(df))

    significance_threshold = 0.00001
    tupple_snarl = df[df['P_value'] < significance_threshold][['Snarl', 'P_value']].itertuples(index=False, name=None)
    write_significative_snarl(tupple_snarl, output_snarl)

    significance_threshold_log = -np.log10(0.00001)
    plt.figure(figsize=(12, 6))
    plt.plot(df['index'], df['-log10(P_value)'], 'o', color='skyblue', alpha=0.7, label='P-value')

    plt.axhline(y=significance_threshold_log, color='red', linestyle='--', label='p=0.00001')
    plt.title('Manhattan Plot of P-values', fontsize=16)
    plt.xlabel('Index', fontsize=14)
    plt.ylabel('-log10(P-value)', fontsize=14)
    plt.legend(loc='upper right')
    plt.savefig(output_file_manh, format='png', dpi=300)

def write_significative_snarl(tupple_snarl, output_snarl) :
    with open(output_snarl, "w") as f :
        for id_snarl, p_value in tupple_snarl :
            f.write(f'{id_snarl}\t{p_value}\n')

if __name__ == "__main__" :

    # file_path = 'output/simulation_1000_binary.tsv'
    # output_file_dist = "output/pgtest_1000_distribution_plot_binary.png"
    # output_file_manh = "output/pgtest_1000_manhattan_plot_binary.png"
    # plot_p_value_distribution_binary(file_path, output_file_dist)
    # plot_manhattan_binary(file_path, output_file_manh)

    # file_path = 'output/simulation_quantitative.tsv'
    # output_file_dist = "output/pgtest_1000_distribution_plot_quantitative.png"
    # output_file_manh = "output/pgtest_1000_manhattan_plot_quantitative.png"
    # plot_p_value_distribution_quantitative(file_path, output_file_dist)
    # plot_manhattan_quantitative(file_path, output_file_manh)

    file_path = 'output/droso_quantitative.tsv'
    output_file_dist = "output/distribution_plot_quantitative.png"
    output_file_manh = "output/manhattan_plot_quantitative.png"
    plot_p_value_distribution_quantitative(file_path, output_file_dist)
    plot_manhattan_quantitative(file_path, output_file_manh)
