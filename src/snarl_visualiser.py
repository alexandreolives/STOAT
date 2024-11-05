import argparse
import subprocess
import json
from cyvcf2 import VCF
import networkx as nx
import os
import pandas as pd


"""
Intuition : 
- Input snarl 
- Make dic : snarl(str) : paths(list) using vg find or list path
- Recalcul group paths
- Output snarl gaf 
"""

def extract_paths_from_subgraph(subgraph_json):
    """
    Parse the subgraph JSON to extract paths between nodes.
    """

    # Create a directed graph
    G = nx.DiGraph()
    for node in subgraph_json.get('nodes', []):
        G.add_node(node['id'])
    for edge in subgraph_json.get('edges', []):
        G.add_edge(edge['from'], edge['to'])
    
    # Find all paths between start and end nodes and return them
    all_paths_dict = {}
    for node in G.nodes():
        # For each node, find paths to all other nodes
        paths = list(nx.all_simple_paths(G, source=node, target=None))
        all_paths_dict[node] = paths
    return all_paths_dict

def get_subgraph(pangenome_graph, start_node, context_size):
    """
    Use `vg find` to get the subgraph as JSON.
    """
    try:
        # Run vg find command and capture the output
        cmd = ['vg', 'find', '-x', pangenome_graph, '-n', str(start_node), '-c', str(context_size)]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        subgraph_json = json.loads(result.stdout)
        return subgraph_json
    except subprocess.CalledProcessError as e:
        print(f"Error executing vg find: {e}")
        return None

def find_sample_paths(snarls: list, vcf_path: str):
    # Dictionary to store which samples pass through which snarl paths
    snarl_matches = {snarl: [] for snarl in snarls}
 
    # Parse VCF file line by line
    for variant in VCF(vcf_path):
        snarl_list = variant.INFO.get('AT', '').split(',')
        genotypes = variant.genotypes  # Extract genotypes once per variant

        # Iterate through the snarl paths
        for snarl_path in snarls:
            # Check if any portion of the snarl path matches the snarl list
            for segment in snarl_list:
                if segment in snarl_path:
                    # Add samples that match this path based on their genotypes
                    for i, genotype in enumerate(genotypes):
                        # Check if the genotype indicates the sample passes through (e.g., not missing)
                        if genotype[0] != -1:  # Replace with an appropriate condition for your use case
                            sample_name = variant.samples[i]
                            if sample_name not in snarl_matches[snarl_path]:
                                snarl_matches[snarl_path].append(sample_name)

    # Return dictionary of snarl paths and matching samples
    return snarl_matches

def parse_group_file(group_file : str):

    df = pd.read_csv(group_file, sep='\t')
    
    # Create dictionaries for group 0 and group 1
    group_0 = {sample: 0 for sample in df[df['GROUP'] == 0]['SAMPLE']}
    group_1 = {sample: 1 for sample in df[df['GROUP'] == 1]['SAMPLE']}
    return group_0, group_1
 
def check_format_pheno_b(file_path : str) -> str :

    if not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")

    with open(file_path, 'r') as file:
        first_line = file.readline().strip()

    header = first_line.split('\t')
    expected_header = ['SAMPLE', 'GROUP']
    if header != expected_header:
        raise argparse.ArgumentTypeError(f"The file must contain the following headers: {expected_header} and be split by tabulation")

    return file_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='List paths through the netgraph of each snarl in a pangenome')
    parser.add_argument('-p', help='the input pangenome .pg file', required=True)
    parser.add_argument('-d', help='the input distance index .dist file', required=True)
    parser.add_argument('-v', help='snarl to decompose', required=False)
    parser.add_argument('-o', help='the output TSV file', type=str, required=True)
    parser.add_argument("-b", "--binary", type=check_format_pheno_b, help="Path to the binary group file (.txt or .tsv)")

    args = parser.parse_args()

    binary_group = parse_group_file(args.binary)

    # Replace these with real start nodes or derive from snarl analysis logic if needed
    start_nodes = [args.v] if args.v else range(1000)  # Example range, adjust as needed
    context_size = 10  # Can be parameterized as needed

    snarl_paths = {}

    for start_node in start_nodes:
        print(f"Processing snarl starting at node {start_node}...")
        subgraph = get_subgraph(args.p, start_node, context_size)
        if subgraph:
            paths = extract_paths_from_subgraph(subgraph)
            snarl_paths[start_node] = paths

    # Write output to TSV
    with open(args.o, 'w') as f:
        for snarl, paths in snarl_paths.items():
            for path in paths.get(snarl, []):
                path_str = '->'.join(map(str, path))
                f.write(f"{snarl}\t{path_str}\n")
    
    print(f"Path extraction completed. Results saved to {args.o}")
