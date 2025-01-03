import argparse
import pandas as pd
from sklearn.metrics import confusion_matrix, precision_score, recall_score, f1_score
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import re

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

        # If group freq diff > threshold, then true label is 0, else 1
        if diff > threshold:
            true_labels.append(0)
        else :
            true_labels.append(1)

        path_list.append(f"{int(row_1['start_node'])}_{int(row_1['next_node'])}")
        list_diff.append(float(diff))
    
    return path_list, true_labels, list_diff

def match_snarl(path_list, true_labels, p_value_file, paths_file, type_):

    p_value_df = pd.read_csv(p_value_file, sep='\t')
    paths_df = pd.read_csv(paths_file, sep='\t')['paths']
    split = p_value_df['SNARL'].str.split('_')

    # To store predicted labels
    predicted_labels_10_2 = []
    predicted_labels_10_5 = []
    predicted_labels_10_8 = []
    cleaned_true_labels = []
    pvalue = []
    num_sample =[]

    for idx, snarl_id in enumerate(path_list):

        start_node, next_node = map(int, snarl_id.split('_'))
        # We want to know if the snarl is in the range of the snarl in the p_value file
        matched_row = p_value_df[(split.str[1].astype(int) <= start_node) & (split.str[0].astype(int) >= next_node) |
                                 (split.str[0].astype(int) <= start_node) & (split.str[1].astype(int) >= next_node)]

        # print("snarl_id  : ", snarl_id)
        # print("matched_row : ", matched_row)

        # Case where the snarl is found 
        if not matched_row.empty:
            indices = matched_row.index
            # print("indices : ", indices)
            split_paths = [paths_df[idx] for idx in indices]
            # print("split_paths : ", split_paths)
            # Check if at least one path in the snarl contains the start node followed by the next node
            for idx_paths, list_path in enumerate(split_paths):
                # print("list_path : ", list_path)
                for path in list_path.split(',') :
                    # print("path : ", path)
                    decomposed_path = [int(num) for num in re.findall(r'\d+', path)]
                    # print("decomposed_path : ", decomposed_path)
                    if start_node in decomposed_path and next_node in decomposed_path:
                        p_value = matched_row.loc[indices[idx_paths]]
                        if type_ == 'binary':  
                            p_value = p_value['P_FISHER']
                        elif type_ == 'quantitative':
                            p_value = p_value['P']
                        else :
                            raise ValueError("type_ must be binary or quantitative")
                        # print("p_value : ", p_value)

                        predicted_labels_10_2.append(0 if p_value < 0.01 else 1)
                        predicted_labels_10_5.append(0 if p_value < 0.00001 else 1)
                        predicted_labels_10_8.append(0 if p_value < 0.00000001 else 1)
                        cleaned_true_labels.append(true_labels[idx])
                        pvalue.append(p_value)
                        num_sample.append(matched_row.loc[indices[idx_paths]]['ALLELE_NUM'])
                        break

    # print("total row : ", len(path_list))
    return predicted_labels_10_2, predicted_labels_10_5, predicted_labels_10_8, cleaned_true_labels, pvalue, num_sample

def conf_mat_maker(p_val, predicted_labels, true_labels, output):
        
    # Calculate confusion matrix for p-value < p_val
    print(f"\nMetrics for p-value < {p_val}:")
    # Inverse because i want the X axis to be the truth labels and Y axis to be the predicted labels
    cm = confusion_matrix(predicted_labels, true_labels)
    print(f"Confusion Matrix for p-value < {p_val}:\n{cm}")
    prec = precision_score(predicted_labels, true_labels)
    recall = recall_score(predicted_labels, true_labels)
    f1 = f1_score(predicted_labels, true_labels)
    print(f"Precision: {prec:.3f}")
    print(f"Recall: {recall:.3f}")
    print(f"F1 Score: {f1:.3f}")

    # Plot confusion matrix for p-value < p_val
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False,
                xticklabels=['Positive', 'Negative'], 
                yticklabels=['Positive', 'Negative'],
                annot_kws={"size": 30})
    plt.xticks(fontsize=16)  
    plt.yticks(fontsize=16)  
    plt.title(f'Confusion Matrix for p-value < {p_val}', fontsize=18)  # Increase title font size
    plt.xlabel('Truth Labels', fontsize=20)  # Increase x-label font size
    plt.ylabel('Predicted Labels', fontsize=20)  # Increase y-label font size
    plt.savefig(output + f'_{p_val}.png', format='png', dpi=300)

def print_confusion_matrix(predicted_labels_10_2, predicted_labels_10_5, predicted_labels_10_8, true_labels, output):
    
    p_val_10_2 = 0.01
    p_val_10_5 = 0.00001
    p_val_10_8 = 0.00000001

    conf_mat_maker(p_val_10_2, predicted_labels_10_2, true_labels, output)
    conf_mat_maker(p_val_10_5, predicted_labels_10_5, true_labels, output)
    conf_mat_maker(p_val_10_8, predicted_labels_10_8, true_labels, output)

def p_value_distribution(test_predicted_labels, cleaned_true_labels, list_diff, p_value, num_sample, output):
    
    false_positive_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 0 and true == 1]

    print("len(false_positive_indices) : " , len(false_positive_indices))

    true_positive_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 0 and true == 0]
    
    print("len(true_positive_indices) : " , len(true_positive_indices))

    false_negative_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 1 and true == 0]

    print("len(false_negative_indices) : " , len(false_negative_indices))

    diff_false_positive = [list_diff[i] for i in false_positive_indices]
    pvalue_false_positive = [p_value[i] for i in false_positive_indices]
    minsample_false_positive = [num_sample[i] for i in false_positive_indices]

    diff_true_positives = [list_diff[i] for i in true_positive_indices]
    pvalue_true_positives = [p_value[i] for i in true_positive_indices]
    minsample_true_positives = [num_sample[i] for i in true_positive_indices]

    diff_false_negative = [list_diff[i] for i in false_negative_indices]
    pvalue_false_negative = [p_value[i] for i in false_negative_indices]
    minsample_false_negative = [num_sample[i] for i in false_negative_indices]

    # Create a DataFrame for easy plotting
    data = {
        'P-Value': pvalue_false_positive + pvalue_true_positives + pvalue_false_negative,
        'Difference': diff_false_positive + diff_true_positives + diff_false_negative,
        'Min Sample': minsample_false_positive + minsample_true_positives + minsample_false_negative,
        'Type': ['False Positives'] * len(pvalue_false_positive) + ['True Positives'] * len(pvalue_true_positives) + ['False Negatives'] * len(pvalue_false_negative)
    }

    df = pd.DataFrame(data)

    # Create the interactive scatter plot
    fig = px.scatter(
        df, 
        x='P-Value', 
        y='Difference', 
        size='Min Sample', 
        color='Type',
        hover_name=df.index,
        title="Distribution of P-Values for False Positives and True Positives",
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
    fig.write_html(f'{output}_pvalue_interactive.html')

    # Extract the corresponding differences for false negatives
    true_negative_diffs = [list_diff[i]*100 for i in false_positive_indices]

    # Plot the distribution of differences for false negatives
    plt.figure(figsize=(10, 6))
    sns.histplot(true_negative_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for False Positives", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig(f'{output}_distribution_false_pos.png', format='png', dpi=300)

def plot_diff_distribution(test_predicted_labels, cleaned_true_labels, list_diff, output, pvalue):

    # True label 0 : difference freq significative
    # True label 1 : No difference freq significative
    # Predict label 0 : pvalue significative
    # Predict label 1 : No pvalue significative

    # ----------------------------- True Positive -----------------------------
    true_positive_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 0 and true == 0]

    true_positive_diffs = [list_diff[i]*100 for i in true_positive_indices]

    plt.figure(figsize=(10, 6))
    sns.histplot(true_positive_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for True Positive", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig(output + pvalue + '_distribution_true_positive.png', format='png', dpi=300)

    # ----------------------------- False Positive -----------------------------
    false_positive_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 0 and true == 1]

    false_positive_diffs = [list_diff[i]*100 for i in false_positive_indices]

    plt.figure(figsize=(10, 6))
    sns.histplot(false_positive_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for False Positive", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig(output + pvalue + '_distribution_false_positive.png', format='png', dpi=300)

    # ----------------------------- False Negative -----------------------------
    false_negative_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 1 and true == 0]

    false_negative_diffs = [list_diff[i]*100 for i in false_negative_indices]

    plt.figure(figsize=(10, 6))
    sns.histplot(false_negative_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for False Negative", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig(output + pvalue + '_distribution_false_negative.png', format='png', dpi=300)

    # ----------------------------- True Negative -----------------------------
    true_negative_indices = [
        i for i, (pred, true) in enumerate(zip(test_predicted_labels, cleaned_true_labels)) 
        if pred == 1 and true == 1]

    true_negative_diffs = [list_diff[i]*100 for i in true_negative_indices]

    plt.figure(figsize=(10, 6))
    sns.histplot(true_negative_diffs, bins=20, kde=True, color='blue')
    plt.title("Distribution of Differences for True Negatives", fontsize=16)
    plt.xlabel("Difference (%)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.grid(False)
    plt.savefig(output + pvalue + '_distribution_true_negatives.png', format='png', dpi=300)
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Parse a file with start_node, next_node, group, and freq columns.")
    parser.add_argument("freq", help="Path to allele frequence file")
    parser.add_argument("p_value", help="Path to p_value file")
    parser.add_argument("paths", help="Path to p_value file")
    parser.add_argument("-t", "--threshold", type=float, required=False, help="Threshold to define the truth label")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--binary", action='store_true', help="binary test")
    group.add_argument("-q", "--quantitative", action='store_true', help="quantitative test")
    args = parser.parse_args()

    if args.binary:
        output = "tests/binary_tests_output/"
        output_diff = f"{output}/binary_test"
        type_ = 'binary'

    elif args.quantitative:
        output = "tests/quantitative_tests_output/"
        output_diff = f"{output}/quantitative_test"
        type_ = 'quantitative'

    # # THRESHOLD_FREQ : Threshold to define the truth label by comparing the frequence difference between both group
    # THRESHOLD_FREQ = 0.4 # Default value is =>20% difference == Truth label 0

    # # Define Truth label from freq file
    # test_path_list, test_true_labels, list_diff = process_file(args.freq, THRESHOLD_FREQ)

    # # Define Truth label from output file and compare with the truth label
    # test_predicted_labels_10_2, test_predicted_labels_10_5, test_predicted_labels_10_8, cleaned_true_labels, _, _ = match_snarl(test_path_list, test_true_labels, args.p_value, args.paths, type_)
        
    # # Plot confusion matrix
    # print_confusion_matrix(test_predicted_labels_10_2, test_predicted_labels_10_5, test_predicted_labels_10_8, cleaned_true_labels, f"{output}/confusion_matrix_{THRESHOLD_FREQ}")

    # THRESHOLD_FREQ = 0.0 : Case where just a difference between both group snarl is considered like Truth label
    THRESHOLD_FREQ = 0.0
    test_path_list, test_true_labels, list_diff = process_file(args.freq, THRESHOLD_FREQ)
    test_predicted_labels_10_2, test_predicted_labels_10_5, test_predicted_labels_10_8, cleaned_true_labels, pvalue, num_sample = match_snarl(test_path_list, test_true_labels, args.p_value, args.paths, type_)
    #print_confusion_matrix(test_predicted_labels_10_2, test_predicted_labels_10_5, test_predicted_labels_10_8, cleaned_true_labels, f"{output}/confusion_matrix_{THRESHOLD_FREQ}")

    # Plot distribution of p-values for false negatives and true positives
    #plot_diff_distribution(test_predicted_labels_10_2, cleaned_true_labels, list_diff, output_diff, "10^-2")
    p_value_distribution(test_predicted_labels_10_2, cleaned_true_labels, list_diff, pvalue, num_sample, output_diff)
    
    """
    python3 tests/verify_truth.py tests/simulation/quantitative_data/pg.snarls.freq.tsv \
    tests/quantitative_tests_output/quantitative_test.assoc.tsv tests/simulation/quantitative_data/snarl_paths.tsv -q
    """

    """
    python3 tests/verify_truth.py tests/simulation/binary_data/pg.snarls.freq.tsv \
    tests/binary_tests_output/binary_test.assoc.tsv tests/simulation/binary_data/snarl_paths.tsv -q
    """