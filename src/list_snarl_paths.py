import bdsg
import argparse
from collections import defaultdict

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

def check_threshold(proportion) :
    proportion = float(proportion)
    if proportion <= 0 or proportion >= 10000 :
        raise ValueError("Proportion value must be >0 and <10000.")

    return proportion

def find_snarl_id(stree, snarl) :
    # create a snarl ID as LEFT_RIGTH bondary nodes
    sstart = stree.get_bound(snarl, False, True)
    sstart = stree.get_node_from_sentinel(sstart)
    send = stree.get_bound(snarl, True, True)
    send = stree.get_node_from_sentinel(send)
    snarl_id = '{}_{}'.format(stree.node_id(sstart),
                            stree.node_id(send))
    
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

def create_snarls(stree, root) :

    # list storing the snarl objects
    snarls = []

    def save_snarl_tree_node(net):
        if stree.is_snarl(net):
            snarls.append(net)

        if not stree.is_node(net) and not stree.is_sentinel(net):
            stree.for_each_child(net, save_snarl_tree_node)
        return (True)
    
    stree.for_each_child(root, save_snarl_tree_node)
    snarls_length = len(snarls)
    print('{} snarls found'.format(snarls_length))

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

def fill_pretty_paths(stree, finished_paths, pretty_paths) :

    for path in finished_paths:
        ppath = Path()
        for net in path:
            if stree.is_sentinel(net):
                net = stree.get_node_from_sentinel(net)
            # Bug enter here with True False
            if stree.is_node(net) or stree.is_trivial_chain(net):
                # if it's a node, add it to the path
                ppath.addNodeHandle(net, stree)

            elif stree.is_chain(net):
                # if it's a chain, we need to write someting like ">Nl>*>Nr"
                nodl = stree.get_bound(net, False, True)
                nodr = stree.get_bound(net, True, False)
                ppath.addNodeHandle(nodl, stree)
                ppath.addNode('*', '>')
                ppath.addNodeHandle(nodr, stree)

        # check if path is mostly traversing nodes in reverse orientation
        if ppath.nreversed() > ppath.size() / 2:
            ppath.flip()
        pretty_paths.append(ppath.print())
    
    return pretty_paths

def write_header_output(output_file) :
    with open(output_file, 'w') as outf:
        outf.write('snarl\tpaths\n')

def write_output(output_file, snarl_id, pretty_paths) :
    with open(output_file, 'a') as outf:
        outf.write('{}\t{}\n'.format(snarl_id, ','.join(pretty_paths)))

def loop_over_snarls_write(stree, snarls, pg, output_file, threshold=50) :
    write_header_output(output_file)

    children = [0]
    def count_children(net):
        children[0] += 1
        return (True)

    # for each snarl, lists paths through the netgraph and write to output TSV
    for idx, snarl in enumerate(snarls):

        if idx > 0 and idx % 10000 == 0 :
            print("idx : ", idx)

        children = [0]
        stree.for_each_child(snarl, count_children)
        if children[0] > threshold :
            print(f"number of children > {threshold}")
            continue

        snarl_id = find_snarl_id(stree, snarl)

        # we'll traverse the netgraph starting at the left boundary
        # init unfinished paths to the first boundary node
        paths = [[stree.get_bound(snarl, False, True)]]
        finished_paths = []
        while len(paths) > 0 :
            path = paths.pop()

            if len(finished_paths) > 10000 :
                print("len of finished_paths > 10000")
                break

            follow_edges(stree, finished_paths, path, paths, pg)

        # prepare path list to output and write each path directly to the file
        pretty_paths = []
        pretty_paths = fill_pretty_paths(stree, finished_paths, pretty_paths)
        write_output(output_file, snarl_id, pretty_paths)

def loop_over_snarls(stree, snarls, pg, threshold=50) :

    snarl_paths = defaultdict(list)

    children = [0]
    def count_children(net):
        children[0] += 1
        return (True)

    # for each snarl, lists paths through the netgraph and write to output TSV
    for snarl in snarls:
        
        if threshold :
            children = [0]
            stree.for_each_child(snarl, count_children)
            if children[0] > threshold :
                #print(f"number of children > {threshold}")
                continue

        snarl_id = find_snarl_id(stree, snarl)
        
        # we'll traverse the netgraph starting at the left boundary
        # init unfinished paths to the first boundary node
        paths = [[stree.get_bound(snarl, False, True)]]
        finished_paths = []
        while len(paths) > 0 :

            path = paths.pop()

            if len(finished_paths) > 10000 :
                #print("len of finished_paths > 10000")
                break

            follow_edges(stree, finished_paths, path, paths, pg)

        # prepare path list to output and write each path directly to the file
        pretty_paths = []
        pretty_paths = fill_pretty_paths(stree, finished_paths, pretty_paths)
        snarl_paths[snarl_id].extend(pretty_paths)

    return snarl_paths

if __name__ == "__main__" :

    parser = argparse.ArgumentParser('List path through the netgraph of each snarl in a pangenome')
    parser.add_argument('-p', help='the input pangenome .pg file', required=True)
    parser.add_argument('-d', help='the input distance index .dist file',required=True)
    parser.add_argument('-t', type=check_threshold, help='Children threshold', required=False)
    parser.add_argument('-o', help='the output TSV file', type=str, required=True)
    args = parser.parse_args()

    stree, pg, root = parse_graph_tree(args.p, args.d)
    snarls = create_snarls(stree, root)

    if args.t :
        loop_over_snarls_write(stree, snarls, pg, args.o, args.t)
    else : 
        loop_over_snarls_write(stree, snarls, pg, args.o)

    # python3 src/list_snarl_paths.py -p ../../snarl_data/fly.pg -d ../../snarl_data/fly.dist -o test_list_snarl.tsv

    # vg find -x ../snarl_data/fly.gbz -r 5176878:5176884 -c 10 | vg view -dp - | dot -Tsvg -o ../snarl_data/subgraph.svg