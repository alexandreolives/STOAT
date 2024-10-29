import argparse
import pandas as pd
from sklearn.metrics import confusion_matrix, precision_score, recall_score, f1_score
import matplotlib.pyplot as plt
import seaborn as sns

# Function to process the frequency file and get result list with differences
def process_file(freq_file, threshold=0.2):
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

    cleaned_true_labels = []
    finded_row = 0

    for idx, snarl_id in enumerate(path_list):

        # if idx == 15 :
        #     exit()
        start_node, next_node = map(int, snarl_id.split('_'))
        split = p_value_df['Snarl'].str.split('_')
        _id_snarl = (split.str[1].astype(int) <= start_node) & (split.str[0].astype(int) >= next_node)
        matched_row = p_value_df[_id_snarl]

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
            p_value_fisher = match['P_value_Fisher']
            predicted_labels_10_2.append(0 if p_value_fisher < 0.01 else 1)
            predicted_labels_10_5.append(0 if p_value_fisher < 0.00001 else 1)
            predicted_labels_10_8.append(0 if p_value_fisher < 0.00000001 else 1)

            finded_row += 1
            cleaned_true_labels.append(true_labels[idx])
    
    print("finded row : ", finded_row)
    print("total row : ", len(path_list))

    return predicted_labels_10_2, predicted_labels_10_5, predicted_labels_10_8, cleaned_true_labels

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

def plot_diff_distribution(test_predicted_labels, cleaned_true_labels, list_diff):
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
    plt.savefig(output + '_distribution_neg_neg.png', format='png', dpi=300)

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
    path_list, true_labels, _ = process_file(args.freq)

    output = "output/truth_confusion_matrix"
    # Match the result_list with the p_value file and write to the output
    predicted_labels_10_2, predicted_labels_10_5, predicted_labels_10_8, cleaned_true_labels = match_snarl(path_list, true_labels, args.p_value, args.paths)
    print_confusion_matrix(predicted_labels_10_2, predicted_labels_10_5, predicted_labels_10_8, cleaned_true_labels, output)

    # special case
    test_path_list, test_true_labels, list_diff = process_file(args.freq, 0.0000000000000001)
    p_val_10_2 = 0.01
    test_predicted_labels_10_2, _, _, cleaned_true_labels = match_snarl(test_path_list, test_true_labels, args.p_value, args.paths)
    conf_mat_maker(p_val_10_2, test_predicted_labels_10_2, cleaned_true_labels, "output/test_truth_confusion_matrix")
    plot_diff_distribution(test_predicted_labels_10_2, cleaned_true_labels, list_diff)

"""
    python3 src/verify_truth.py ../../snarl_data/simulation_1000vars_100samps/pg.snarls.freq.tsv output/simulation_1000_binary.tsv
"""
