import bdsg
import argparse
import re
from collections import defaultdict
import time 

# class to help make paths from BDSG objects
# and deal with orientation, flipping, etc
class Path:
    def __init__(self):
        self.nodes = []
        self.orients = []

    def addNode(self, node, orient):
        # here we know the actual node id and orientation
        self.nodes.append(node)
        self.orients.append(orient)

    def addNodeHandle(self, node_h, stree):
        # we have a BDSG handle and need to extract node id and orientation
        # couldn't find a better way than parsing the
        # string representation of the node...
        node_s = stree.net_handle_as_string(node_h)
        # trivial chain
        if stree.is_trivial_chain(node_h):
            node_s = node_s.replace(' pretending to be a chain', '')
            node_s = node_s.replace(' in a simple snarl', '')

        # parse node info
        node_s = node_s.replace('node ', '')
        node_o = '>'
        if 'rev' in node_s:
            node_o = '<'
        node_s = node_s.replace('rev', '').replace('fd', '')
        # add node to path
        self.nodes.append(node_s)
        self.orients.append(node_o)

    def print(self):
        # write the string representation of the path
        # e.g. ">I>J<K>L" or ">I>J>*>M>N"
        out_path = []
        for ii in range(len(self.nodes)):
            out_path.append(self.orients[ii] + str(self.nodes[ii]))
        return (''.join(out_path))

    def flip(self):
        # some paths are traversing the nodes entirely or mostly in reverse
        # for convenience we might can to flip them
        self.nodes.reverse()
        self.orients.reverse()
        for ii in range(len(self.orients)):
            if self.nodes[ii] == '*':
                # don't flip the orientation of the * because
                # it should be ">" by construction
                continue
            if self.orients[ii] == '>':
                self.orients[ii] = '<'
            else:
                self.orients[ii] = '>'

    def size(self):
        return (len(self.nodes))

    def nreversed(self):
        # counts how many nodes are traversed in reverse
        return (sum(['<' == orient for orient in self.orients]))

def split_paths(path) :
    return re.findall(r'\d+', path)

def length_node(pg, node_id) :
    return pg.get_length(node_id)

# def calcul_type_variant(stree, type_variants) :
#     """ 
#     Calcul the type variant of a tested snarl
#     If snarl are only node of length 1 => SNP 
#     Else => COMPLEX (unknow variant type)
#     """
#     list_type_variant = []
#     for length in type_variants :
#         list_node = split_paths(path)
#         if len(list_node) == 3 : # Case simple path len 3
#             length_middle_node = length_node(stree, list_node[1])
#             list_type_variant.append(str(length_middle_node))

#         if len(list_node) > 3 : # Case snarl in snarl / Indel
#             list_type_variant.append("COMPLEX")

#         else : # Deletion
#             list_type_variant.append("0")

#     return list_type_variant

def check_threshold(proportion) :
    proportion = float(proportion)
    if proportion <= 0 :
        raise argparse.ArgumentTypeError("Proportion value must be >0.")

    return proportion

def find_snarl_id(stree, snarl) :
    # create a snarl ID as LEFT_RIGTH bondary nodes
    sstart = stree.get_bound(snarl, False, True)
    sstart = stree.get_node_from_sentinel(sstart)
    send = stree.get_bound(snarl, True, True)
    send = stree.get_node_from_sentinel(send)
    snarl_id = '{}_{}'.format(stree.node_id(send), stree.node_id(sstart))

    return snarl_id

def follow_edges(stree, finished_paths, path, paths, pg) :
    def add_to_path(next_child) :

        if stree.is_sentinel(next_child):
            # If this is the bound of the snarl then we're done
            # Because we only traverse in the netgraph, it can only be the
            # bound of the parent snarl
            finished_paths.append([])
            for net in path:
                finished_paths[-1].append(net)
            finished_paths[-1].append(next_child)
        else :
            for i in path : 
                # Case where we find a loop 
                if stree.net_handle_as_string(i) == stree.net_handle_as_string(next_child) :
                    return False
                
            paths.append([])
            for net in path:
                paths[-1].append(net) #paths[-1].append(net)
            paths[-1].append(next_child) #paths[-1].append(next_child)
        return True

    # from the last thing in the path
    stree.follow_net_edges(path[-1], pg, False, add_to_path)

def save_snarls(stree, root) :

    # list storing the snarl objects
    snarls = []

    def save_snarl_tree_node(net):
        if stree.is_snarl(net):
            snarls.append(net)

        if not stree.is_node(net) and not stree.is_sentinel(net):
            stree.for_each_child(net, save_snarl_tree_node)
        return (True)
    
    stree.for_each_child(root, save_snarl_tree_node)
    return snarls

def parse_graph_tree(pg_file, dist_file) :

    # load graph and snarl tree
    pg = bdsg.bdsg.PackedGraph()
    pg.deserialize(pg_file)
    stree = bdsg.bdsg.SnarlDistanceIndex()
    stree.deserialize(dist_file)

    # list all snarls in pangenome
    # init with the child (only one ideally) of the root
    root = stree.get_root()
    return stree, pg, root

def fill_pretty_paths(stree, pg, finished_paths) :
    pretty_paths = []
    length_net_paths = []

    for path in finished_paths:
        ppath = Path()
        length_net = []
        for net in path:
            if stree.is_sentinel(net):
                net = stree.get_node_from_sentinel(net)

            if stree.is_node(net) or stree.is_trivial_chain(net):
                # if it's a node, add it to the path
                ppath.addNodeHandle(net, stree)
                if stree.is_node(net) :
                    length_net.append(str(stree.node_length(net)))

                else :
                    stn_start = stree.get_bound(net, False, True)
                    stree.node_id(stn_start)
                    net_trivial_chain = stree.get_handle(stn_start, pg)
                    length_net.append(str(stree.node_length(net_trivial_chain)))

            elif stree.is_chain(net):
                # if it's a chain, we need to write someting like ">Nl>*>Nr"
                nodl = stree.get_bound(net, False, True)
                nodr = stree.get_bound(net, True, False)
                ppath.addNodeHandle(nodl, stree)
                ppath.addNode('*', '>')
                ppath.addNodeHandle(nodr, stree)
                length_net.append("-1")

        # check if path is mostly traversing nodes in reverse orientation
        if ppath.nreversed() > ppath.size() / 2:
            ppath.flip()
        pretty_paths.append(ppath.print())
        length_net_paths.extend(length_net)

    return pretty_paths, length_net_paths

def write_header_output(output_file) :
    with open(output_file, 'w') as outf:
        outf.write('snarl\tpaths\ttype\n')

def write_output(output_file, snarl_id, pretty_paths, type_variants) :
    with open(output_file, 'a') as outf:
        outf.write('{}\t{}\t{}\n'.format(snarl_id, ','.join(pretty_paths), ','.join(type_variants)))

def write_header_output_not_analyse(output_file) :
    with open(output_file, 'w') as outf:
        outf.write('snarl\treason\n')

def write_output_not_analyse(output_file, snarl_id, reason) :
    with open(output_file, 'a') as outf:
        outf.write('{}\t{}\n'.format(snarl_id, reason))

def loop_over_snarls_write(stree, snarls, pg, output_file, output_snarl_not_analyse, time_threshold=10) :

    write_header_output(output_file)
    write_header_output_not_analyse(output_snarl_not_analyse)

    # for each snarl, lists paths through the netgraph and write to output TSV
    for snarl in snarls:

        snarl_time = time.time()
        snarl_id = find_snarl_id(stree, snarl)

        # we'll traverse the netgraph starting at the left boundary
        # init unfinished paths to the first boundary node
        paths = [[stree.get_bound(snarl, False, True)]]
        finished_paths = []
        while len(paths) > 0 :
            path = paths.pop()

            if snarl_time - time.time() > time_threshold :
                write_output_not_analyse(output_snarl_not_analyse, snarl_id, "time_calculation_out")
                break

            follow_edges(stree, finished_paths, path, paths, pg)

        # prepare path list to output and write each path directly to the file
        pretty_paths, type_variants = fill_pretty_paths(stree, pg, finished_paths)
        print("pretty_paths : ", pretty_paths)
        print("type_variants : ", type_variants)
        # type_variants_2 = calcul_type_variant(pg, type_variants)
        # print("type_variants_2 : ", type_variants_2)
        exit()
        write_output(output_file, snarl_id, pretty_paths, type_variants)

def loop_over_snarls(stree, snarls, pg, output_snarl_not_analyse, threshold=50) :

    write_header_output_not_analyse(output_snarl_not_analyse)

    snarl_paths = defaultdict(list)
    snarl_number_analysis = 0

    children = [0]
    def count_children(net):
        children[0] += 1
        return (True)

    # for each snarl, lists paths through the netgraph and write to output TSV
    for snarl in snarls:
        
        children = [0]
        snarl_id = find_snarl_id(stree, snarl)

        stree.for_each_child(snarl, count_children)
        if children[0] > threshold :
            write_output_not_analyse(output_snarl_not_analyse, snarl_id, "too_many_children")
            continue
        
        # we'll traverse the netgraph starting at the left boundary
        # init unfinished paths to the first boundary node
        paths = [[stree.get_bound(snarl, False, True)]]
        finished_paths = []
        while len(paths) > 0 :

            path = paths.pop()

            if len(finished_paths) > 10000 :
                write_output_not_analyse(output_snarl_not_analyse, snarl_id, "number_of_paths_to_hight")
                break

            follow_edges(stree, finished_paths, path, paths, pg)

        # prepare path list to output and write each path directly to the file
        pretty_paths = fill_pretty_paths(stree, pg, finished_paths)
        snarl_paths[snarl_id].extend(pretty_paths)
        snarl_number_analysis += len(pretty_paths)

    return snarl_paths, snarl_number_analysis

if __name__ == "__main__" :

    parser = argparse.ArgumentParser('List path through the netgraph of each snarl in a pangenome')
    parser.add_argument('-p', help='the input pangenome .pg file', required=True)
    parser.add_argument('-d', help='the input distance index .dist file',required=True)
    parser.add_argument('-t', type=check_threshold, help='time calculation threshold', required=False)
    parser.add_argument('-o', help='the output TSV file', type=str, required=True)
    args = parser.parse_args()

    stree, pg, root = parse_graph_tree(args.p, args.d)
    snarls = save_snarls(stree, root)
    print(f"Total of snarls found : {len(snarls)}")
    print("Saving snarl path decomposition...")

    output_snarl_not_analyse = "snarl_not_analyse.tsv"

    threshold = args.t if args.t else 10
    loop_over_snarls_write(stree, snarls, pg, args.o, output_snarl_not_analyse, threshold)

    # python3 src/list_snarl_paths.py -p /home/mbagarre/Bureau/droso_data/fly/fly.pg -d /home/mbagarre/Bureau/droso_data/fly/fly.dist -o output/test/test_list_snarl.tsv
    # vg find -x ../snarl_data/fly.gbz -r 5176878:5176884 -c 10 | vg view -dp - | dot -Tsvg -o ../snarl_data/subgraph.svg
    