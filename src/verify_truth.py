import argparse
import pandas as pd
from sklearn.metrics import confusion_matrix, accuracy_score, recall_score, f1_score
import matplotlib.pyplot as plt
import seaborn as sns

# Function to process the frequency file and get result list with differences
def process_file(file_path, threshold=0.4):
    df = pd.read_csv(file_path, sep='\t')
    path_list = []
    true_labels = []

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

    return path_list, true_labels

def match_snarl(path_list, true_labels, p_value_file):
    p_value_df = pd.read_csv(p_value_file, sep='\t')
    
    # To store predicted labels
    predicted_labels_10_2 = []
    predicted_labels_10_5 = []
    cleaned_true_labels = []
    
    finded_row = 0

    for idx, snarl_id in enumerate(path_list):
        start_node, next_node = map(int, snarl_id.split('_'))

        matched_row = p_value_df[(p_value_df['Snarl'].str.split('_').str[1].astype(int) == start_node) & 
                                    (p_value_df['Snarl'].str.split('_').str[0].astype(int) >= next_node)]

        if not matched_row.empty:
            match = matched_row.iloc[0]
            p_value_fisher = match['P_value_Fisher']

            predicted_labels_10_2.append(0 if p_value_fisher < 0.01 else 1)
            predicted_labels_10_5.append(0 if p_value_fisher < 0.00001 else 1)

            finded_row += 1
            cleaned_true_labels.append(true_labels[idx])
    
    print("finded row : ", finded_row)
    print("total row : ", len(path_list))

    return predicted_labels_10_2, predicted_labels_10_5, cleaned_true_labels

def print_confusion_matrix(predicted_labels_10_2, predicted_labels_10_5, true_labels, output):
    # Calculate confusion matrix for p-value < 0.01
    print("\nMetrics for p-value < 0.01:")
    cm_10_2 = confusion_matrix(true_labels, predicted_labels_10_2)
    print(f"Confusion Matrix for p-value < 0.01:\n{cm_10_2}")
    acc_10_2 = accuracy_score(true_labels, predicted_labels_10_2)
    recall_10_2 = recall_score(true_labels, predicted_labels_10_2)
    f1_10_2 = f1_score(true_labels, predicted_labels_10_2)
    print(f"Accuracy: {acc_10_2:.2f}")
    print(f"Recall: {recall_10_2:.2f}")
    print(f"F1 Score: {f1_10_2:.2f}")

    # Plot confusion matrix for p-value < 0.01
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm_10_2, annot=True, fmt='d', cmap='Blues', cbar=False,
                xticklabels=['Predicted Pos', 'Predicted Neg'], 
                yticklabels=['True Pos', 'True Neg'],
                annot_kws={"size": 30})
    plt.xticks(fontsize=16)  
    plt.yticks(fontsize=16)  
    plt.title('Confusion Matrix for p-value < 0.01', fontsize=18)  # Increase title font size
    plt.xlabel('Predicted Labels', fontsize=20)  # Increase x-label font size
    plt.ylabel('True Labels', fontsize=20)  # Increase y-label font size
    plt.savefig(output + '_0.01.png', format='png', dpi=300)

    # Calculate confusion matrix for p-value < 0.00001
    print("\nMetrics for p-value < 0.00001:")
    cm_10_5 = confusion_matrix(true_labels, predicted_labels_10_5)
    print(f"Confusion Matrix for p-value < 0.00001:\n{cm_10_5}")
    acc_10_5 = accuracy_score(true_labels, predicted_labels_10_5)
    recall_10_5 = recall_score(true_labels, predicted_labels_10_5)
    f1_10_5 = f1_score(true_labels, predicted_labels_10_5)
    print(f"Accuracy: {acc_10_5:.2f}")
    print(f"Recall: {recall_10_5:.2f}")
    print(f"F1 Score: {f1_10_5:.2f}")

    # Plot confusion matrix for p-value < 0.00001
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm_10_5, annot=True, fmt='d', cmap='Blues', cbar=False,
                xticklabels=['Predicted Pos', 'Predicted Neg'], 
                yticklabels=['True Pos', 'True Neg'],
                annot_kws={"size": 30})
    plt.xticks(fontsize=16)  
    plt.yticks(fontsize=16)  
    plt.title('Confusion Matrix for p-value < 0.00001', fontsize=18)
    plt.xlabel('Predicted Labels', fontsize=20)  
    plt.ylabel('True Labels', fontsize=20)  
    plt.savefig(output + '_0.00001.png', format='png', dpi=300)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse a file with start_node, next_node, group, and freq columns.")
    parser.add_argument("freq", help="Path to allele frequence file")
    parser.add_argument("p_value", help="Path to p_value file")

    args = parser.parse_args()
    path_list, true_labels = process_file(args.freq)

    output = "output/truth_confusion_matrix"
    # Match the result_list with the p_value file and write to the output
    predicted_labels_10_2, predicted_labels_10_5, cleaned_true_labels = match_snarl(path_list, true_labels, args.p_value)
    print_confusion_matrix(predicted_labels_10_2, predicted_labels_10_5, cleaned_true_labels, output)

"""
    python3 src/verify_truth.py ../../snarl_data/simulation_1000vars_100samps/pg.snarls.freq.tsv output/simulation_1000_binary.tsv

    finded row :  1532
    total row :  3064

    Metrics for p-value < 0.01:
    Confusion Matrix for p-value < 0.01:
    [[  78    8]
    [  36 1410]]
    Accuracy: 0.9712793733681462
    Recall: 0.975103734439834
    F1 Score: 0.9846368715083799

    Metrics for p-value < 0.00001:
    Confusion Matrix for p-value < 0.00001:
    [[  61   25]
    [   5 1441]]
    Accuracy: 0.9804177545691906
    Recall: 0.9965421853388658
    F1 Score: 0.9896978021978022
"""
