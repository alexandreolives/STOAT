import argparse
import list_snarl_paths
import snarl_vcf_parser

parser = argparse.ArgumentParser('List path through the netgraph of each snarl in a pangenome')
parser.add_argument('-p', help='the input pangenome .pg file', required=True)
parser.add_argument('-d', help='the input distance index .dist file',required=True)
parser.add_argument("-t", "--threshold", type=list_snarl_paths.check_threshold, help='Children threshold', required=False)

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-b", "--binary", type=snarl_vcf_parser.check_format_group_snarl, help="Path to the binary group file (.txt or .tsv)")
group.add_argument("-q", "--quantitative", type=snarl_vcf_parser.check_format_pheno, help="Path to the quantitative phenotype file (.txt or .tsv)")
parser.add_argument("-o", "--output", type=str, required=False, help="Path to the output file")
args = parser.parse_args()

stree, pg, root = list_snarl_paths.parse_graph_tree(args.p, args.d)
snarls = list_snarl_paths.create_snarls(stree, root)

snarl_paths = list_snarl_paths.loop_over_snarls(stree, snarls, pg, args.t)
vcf_object = snarl_vcf_parser.SnarlProcessor(args.vcf_path)
vcf_object.fill_matrix()

if args.binary:
    binary_group = snarl_vcf_parser.parse_group_file(args.binary)
    if args.output :
        vcf_object.binary_table(snarl_paths, binary_group, args.output)
    else :
        vcf_object.binary_table(snarl_paths, binary_group)

if args.quantitative:
    quantitative = snarl_vcf_parser.parse_pheno_file(args.quantitative)
    if args.output :
        vcf_object.quantitative_table(snarl_paths, quantitative, args.output)
    else :
        vcf_object.quantitative_table(snarl_paths, quantitative)
