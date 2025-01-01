import bdsg # type: ignore
import argparse
from collections import Counter

def count_children(net):
    list_unordered_children[0] += 1
    return (True)

def save_snarl_tree_node(net):
    def count_children_in(net):
        list_unordered_children[-1] += 1
        return (True)

    if stree.is_snarl(net):
        list_unordered_children.append(0)
        stree.for_each_child(net, count_children_in)
        snarls.append(net)

    if not stree.is_node(net) and not stree.is_sentinel(net):
        stree.for_each_child(net, save_snarl_tree_node)
    return (True)

def check_threshold(proportion) :
    proportion = float(proportion)
    if proportion <= 0.0 or proportion >= 100.0 :
        raise ValueError("Proportion value must be >0 and <100.")

    return proportion

def calcul(list_unordered_children) :
    # Count the frequency of each number
    frequency = Counter(list_unordered_children)
    sorted_frequency = sorted(frequency.items(), key=lambda x: x[0], reverse=True)
    ordered_childrens, effs = zip(*sorted_frequency)
    childrens, effs = list(ordered_childrens), list(effs)

    children_multiply_by_eff = [children * eff for children, eff in zip(childrens, effs)]
    sum_eff = sum(children_multiply_by_eff)
    sum_cur = 0
    cumulated_proportion = []

    for eff, children in zip(effs, childrens) :
        number_of_children = eff*children
        proportion = (number_of_children + sum_cur)/sum_eff
        cumulated_proportion.append((1-proportion)*100)
        sum_cur += number_of_children
    
    return childrens, effs, cumulated_proportion

def children_cut_off(childrens, effs, cumulated_proportion, proportion_threshold) -> int :
    soon_dead_child = 0
    i = 0
    while cumulated_proportion[i] >= proportion_threshold :
        number_of_children = effs[i]*childrens[i]
        soon_dead_child += number_of_children
        i += 1

    return soon_dead_child

if __name__ == "__main__" :

    parser = argparse.ArgumentParser('Identify the number of children per snarl')
    parser.add_argument('-p', help='the input pangenome .pg file', required=True)
    parser.add_argument('-d', help='the input distance index .dist file',required=True)
    parser.add_argument("-c", "--cut", type=check_threshold, help="print the number of children cut off after applyed the threshold proportion", required=False)
    parser.add_argument("-o", "--output", help='the output TSV file', required=False)
    args = parser.parse_args()

    # load graph and snarl tree
    pg = bdsg.bdsg.PackedGraph()
    pg.deserialize(args.p)
    stree = bdsg.bdsg.SnarlDistanceIndex()
    stree.deserialize(args.d)

    # list all snarls in pangenome
    # init with the child (only one ideally) of the root
    root = stree.get_root()
    snarls = []
    list_unordered_children = [0]

    stree.for_each_child(root, count_children)
    stree.for_each_child(root, save_snarl_tree_node)
    print('{} children found'.format(len(list_unordered_children)))

    childrens, effs, cumulated_proportion = calcul(list_unordered_children)

    if args.output :
        with open(args.output, 'wt') as outf:
            outf.write('Number_children\tEffectif\tCumulated_proportion(%)\n')
            for child, eff, cumul_eff in zip(childrens, effs, cumulated_proportion) :
                outf.write('{}\t{}\t{}\n'.format(child, eff, cumul_eff))

    if args.cut :
        number_dead_children = children_cut_off(childrens, effs, cumulated_proportion, args.cut)
        print('{} children cut off after applyed the threshold proportion'.format(number_dead_children))

    # python3 src/children_distribution.py -p ../../snarl_data/fly.pg -d ../../snarl_data/fly.dist -o output/children_distribution.tsv
