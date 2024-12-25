import argparse
import pandas as pd
from sklearn.metrics import confusion_matrix, precision_score, recall_score, f1_score
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px

# Function to process the frequency file and get result list with differences
def process_file(freq_file, threshold=0.01):
    df = pd.read_csv(freq_file, sep='\t')
    path_list = []
    true_labels = []
    list_diff = []

    # Iterate in pairs and calculate the differences
    for i in range(0, len(df) - 1, 2):
        row_1 = df.iloc[i]
        row_2 = df.iloc[i + 1]
        diff = abs(row_1['freq'] - row_2['freq'])

        if diff >= threshold:
            true_labels.append(0)
        else :
            true_labels.append(1)

        path_list.append(f"{int(row_1['next_node'])}_{int(row_1['start_node'])}")
        list_diff.append(float(diff))

    return path_list, true_labels, list_diff

def match_snarl(path_list, true_labels, p_value_file, paths_file):
    p_value_df = pd.read_csv(p_value_file, sep='\t')
    paths_df = pd.read_csv(paths_file, sep='\t')

    # To store predicted labels
    predicted_labels_10_2 = []
    predicted_labels_10_5 = []
    predicted_labels_10_8 = []
    list_pvalue = []
    list_min_sample = []

    cleaned_true_labels = []
    finded_row = 0

    for idx, snarl_id in enumerate(path_list):

        start_node, next_node = map(int, snarl_id.split('_'))
        split = p_value_df['SNARL'].str.split('_')
        _id_snarl = (split.str[1].astype(int) <= start_node) & (split.str[0].astype(int) >= next_node)
        matched_row = p_value_df[_id_snarl]
        min_row = p_value_df['Min_sample']

        if not matched_row.empty :
            if len(matched_row) > 1 :
                split_paths = paths_df['paths'][_id_snarl]
                indices = split_paths.index
                for idx_paths, list_paths in enumerate(split_paths) :
                    for path in list_paths.split(',') :
                        if f"{next_node}>{start_node}" in path :
                            match = matched_row.loc[indices[idx_paths]]

            else :
                match = matched_row.iloc[0]

            p_value_fisher = match['P_Fisher']
            min_row = match['Min_sample']
            list_pvalue.append(p_value_fisher)
            list_min_sample.append(min_row)
            predicted_labels_10_2.append(0 if p_value_fisher < 0.01 else 1)
            predicted_labels_10_5.append(0 if p_value_fisher < 0.00001 else 1)
            predicted_labels_10_8.append(0 if p_value_fisher < 0.00000001 else 1)

            finded_row += 1
            cleaned_true_labels.append(true_labels[idx])
    
    print("finded row : ", finded_row)
    print("total row : ", len(path_list))

    return predicted_labels_10_2, predicted_labels_10_5, predicted_labels_10_8, cleaned_true_labels, list_pvalue, list_min_sample

def conf_mat_maker(p_val, predicted_labels, true_labels, output) :
    
    # Calculate confusion matrix for p-value < p_val
    print(f"\nMetrics for p-value < {p_val}:")
    cm = confusion_matrix(true_labels, predicted_labels)
    print(f"Confusion Matrix for p-value < {p_val}:\n{cm}")
    prec = precision_score(true_labels, predicted_labels)
    recall = recall_score(true_labels, predicted_labels)
    f1= f1_score(true_labels, predicted_labels)
    print(f"Precision: {prec:.3f}")
    print(f"Recall: {recall:.3f}")
    print(f"F1 Score: {f1:.3f}")

    # Plot confusion matrix for p-value < p_val
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False,
                xticklabels=['Predicted Pos', 'Predicted Neg'], 
                yticklabels=['True Pos', 'True Neg'],
                annot_kws={"size": 30})
    plt.xticks(fontsize=16)  
    plt.yticks(fontsize=16)  
    plt.title(f'Confusion Matrix for p-value < {p_val}', fontsize=18)  # Increase title font size
    plt.xlabel('Predicted Labels', fontsize=20)  # Increase x-label font size
    plt.ylabel('True Labels', fontsize=20)  # Increase y-label font size
    plt.savefig(output + f'_{p_val}.png', format='png', dpi=300)

def print_confusion_matrix(predicted_labels_10_2, predicted_labels_10_5, predicted_labels_10_8, true_labels, output):
    
    p_val_10_2 = 0.01
    p_val_10_5 = 0.00001
    p_val_10_8 = 0.00000001

    conf_mat_maker(p_val_10_2, predicted_labels_10_2, true_labels, output)
    conf_mat_maker(p_val_10_5, predicted_labels_10_5, true_labels, output)
    conf_mat_maker(p_val_10_8, predicted_labels_10_8, true_labels, output)

def p_value_distribution(test_predicted_labels, cleaned_true_labels, list_diff, p_value, min_sample):
    # Identify indices for false negatives and true positives
    false_negatives_indices = [
        i for i, (pred, true, diff) in enumerate(zip(test_predicted_labels, cleaned_true_labels, list_diff))
        if (pred == 1 and true == 0) and diff > 0
    ]

    true_positive_indices = [
        i for i, (pred, true, diff) in enumerate(zip(test_predicted_labels, cleaned_true_labels, list_diff))
        if (pred == 0 and true == 0) and diff > 0
    ]

    # Extract data for false negatives
    diff_false_negatives = [list_diff[i] for i in false_negatives_indices]
    pvalue_false_negatives = [p_value[i] for i in false_negatives_indices]
    minsample_false_negatives = [min_sample[i] for i in false_negatives_indices]

    # Extract data for true positives
    diff_true_positives = [list_diff[i] for i in true_positive_indices]
    pvalue_true_positives = [p_value[i] for i in true_positive_indices]
    minsample_true_positives = [min_sample[i] for i in true_positive_indices]

    # Create a DataFrame for easy plotting
    data = {
        'P-Value': pvalue_false_negatives + pvalue_true_positives,
        'Difference': diff_false_negatives + diff_true_positives,
        'Min Sample': minsample_false_negatives + minsample_true_positives,
        'Type': ['False Negatives'] * len(pvalue_false_negatives) + ['True Positives'] * len(pvalue_true_positives)
    }

    df = pd.DataFrame(data)

    # Create the interactive scatter plot
    fig = px.scatter(
        df, 
        x='P-Value', 
        y='Difference', 
        size='Min Sample', 
        color='Type',
        hover_name=df.index,  # Show the index on hover
        title="Distribution of P-Values for False Negatives and True Positives",
        labels={"P-Value": "P-Value", "Difference": "Simulated Effect (Difference in Probabilities)"},
        size_max=20
    )

    fig.update_layout(
        xaxis_title="P-Value",
        yaxis_title="Simulated Effect (Difference in Probabilities)",
        legend_title="Type",
        template="plotly_white"
    )

    # Show the interactive plot
    fig.show()

    # Optionally, save the plot as an HTML file
    fig.write_html('output/test_pvalue_interactive.html')

    # Extract the corresponding differences for false negatives
    true_negative_diffs = [list_diff[i]*100 for i in false_negatives_indices]

    # Plot the distribution of differences for false negatives
    plt.figure(figsize=(10, 6))
    sns.histplot(true_negative_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for False Negatives", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig('output/test_distribution_false_neg.png', format='png', dpi=300)

def plot_diff_distribution(test_predicted_labels, cleaned_true_labels, list_diff, output):
    # Identify false negatives: True label is 1 (no difference), but predicted label is 0 (predicted a difference)
    true_negatives_indices = [i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) if pred == 1 and true == 0]

    # Extract the corresponding differences for false negatives
    true_negative_diffs = [list_diff[i]*100 for i in true_negatives_indices]

    # Plot the distribution of differences for false negatives
    plt.figure(figsize=(10, 6))
    sns.histplot(true_negative_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for True Negatives", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig(output + '_distribution_true_neg.png', format='png', dpi=300)

    true_true_indices = [i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) if pred == 0 and true == 0]
    true_true_diffs = [list_diff[i]*100 for i in true_true_indices]

    # Plot the distribution of differences for false negatives
    plt.figure(figsize=(10, 6))
    sns.histplot(true_true_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for True True", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig(output + '_distribution_true_true.png', format='png', dpi=300)
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Parse a file with start_node, next_node, group, and freq columns.")
    parser.add_argument("freq", help="Path to allele frequence file")
    parser.add_argument("p_value", help="Path to p_value file")
    parser.add_argument("paths", help="Path to p_value file")
    args = parser.parse_args()

    # Case where all element is a Truth label 
    THRESHOLD = 0.0

    # Define Truth label from freq file
    test_path_list, test_true_labels, list_diff = process_file(args.freq, THRESHOLD)

    # Define Truth label from output file and compare with the truth label
    test_predicted_labels_10_2, test_predicted_labels_10_5, test_predicted_labels_10_8, cleaned_true_labels, list_pvalue, list_min_sample = match_snarl(test_path_list, test_true_labels, args.p_value, args.paths)
        
    # Match the result_list with the p_value file and write to the output
    print_confusion_matrix(test_predicted_labels_10_2, test_predicted_labels_10_5, test_predicted_labels_10_8, cleaned_true_labels, "tests/binary_tests_output/truth_confusion_matrix")

    p_val_10_2 = 0.01
    conf_mat_maker(p_val_10_2, test_predicted_labels_10_2, cleaned_true_labels, "tests/binary_tests_output/binary_test_truth_confusion_matrix")
    plot_diff_distribution(test_predicted_labels_10_2, cleaned_true_labels, list_diff, "tests/binary_tests_output/binary_test")

    p_val_10_5 = 0.00001
    conf_mat_maker(p_val_10_5, test_predicted_labels_10_5, cleaned_true_labels, "tests/binary_tests_output/binary_test_truth_confusion_matrix")
    plot_diff_distribution(test_predicted_labels_10_5, cleaned_true_labels, list_diff, "tests/binary_tests_output/binary_test")

    p_val_10_8 = 0.000000001
    conf_mat_maker(p_val_10_8, test_predicted_labels_10_8, cleaned_true_labels, "tests/binary_tests_output/binary_test_truth_confusion_matrix")
    plot_diff_distribution(test_predicted_labels_10_8, cleaned_true_labels, list_diff, "tests/binary_tests_output/binary_test")

    """
    python3 src/verify_truth.py tests/simulation/pg.snarls.freq.tsv \
    tests/binary_tests_output/binary_test.assoc.tsv tests/simulation/snarl_paths.tsv
    """
