import bdsg
import argparse

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
        # trivial chain?
        if stree.is_trivial_chain(node_h):
            node_s = node_s.replace(' pretending to be a chain '
                                    'in a simple snarl', '')
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

parser = argparse.ArgumentParser('List path through the netgraph of each snarl in a pangenome')
parser.add_argument('-p', help='the input pangenome .pg file', required=True)
parser.add_argument('-d', help='the input distance index .dist file',required=True)
parser.add_argument('-o', help='the output TSV file', required=True)
args = parser.parse_args()

# load graph and snarl tree
pg = bdsg.bdsg.PackedGraph()
pg.deserialize(args.p)
stree = bdsg.bdsg.SnarlDistanceIndex()
stree.deserialize(args.d)

# list all snarls in pangenome
# init with the child (only one ideally) of the root
root = stree.get_root()
# list storing the snarl objects
snarls = []

def save_snarl_tree_node(net):
    if stree.is_snarl(net):
        snarls.append(net)
    if not stree.is_node(net) and not stree.is_sentinel(net):
        stree.for_each_child(net, save_snarl_tree_node)
    return (True)

stree.for_each_child(root, save_snarl_tree_node)
print('{} snarls found'.format(len(snarls)))

# output file will be a TSV file with two columns
outf = open(args.o, 'wt')
outf.write('snarl\tpaths\n')
# count paths, just for log purpose
npaths = 0

# for each snarl, lists paths through the netgraph and write to output TSV
for snarl in snarls:
    # create a snarl ID as LEFT_RIGTH bondary nodes
    sstart = stree.get_bound(snarl, False, True)
    sstart = stree.get_node_from_sentinel(sstart)
    send = stree.get_bound(snarl, True, True)
    send = stree.get_node_from_sentinel(send)
    snarl_id = '{}_{}'.format(stree.node_id(sstart),
                              stree.node_id(send))
    # we'll traverse the netgraph starting at the left boundary
    # init unfinished paths to the first boundary node
    paths = [[stree.get_bound(snarl, False, True)]]
    finished_paths = []
    # while there are still unfinished paths, pop one path, follow net edges
    # to get the next step from the last node in the path, and then either
    # remember the completed path or remember to continue traversing
    while len(paths) > 0:
        path = paths.pop()

        # helper function to add the next child to the path we're building and
        # either add it to the list of completed paths or the list of paths
        # to continue building
        def add_to_path(next_child):
            if stree.is_sentinel(next_child):
                # If this is the bound of the snarl then we're done
                # Because we only traverse in the netgraph, it can only be the
                # bound of the parent snarl
                # Explicitly make a deep copy of the path because
                # idk how to do it in python (<- Xian)
                finished_paths.append([])
                for net in path:
                    finished_paths[-1].append(net)
                finished_paths[-1].append(next_child)
            else:
                # If we reached sibling child of the snarl,
                # then continue the traversal
                # Do another copy into the other list of paths
                paths.append([])
                for net in path:
                    paths[-1].append(net)
                paths[-1].append(next_child)
            return True
        # run add_to_path for everything one step out
        # from the last thing in the path
        stree.follow_net_edges(path[-1], pg, False, lambda n: add_to_path(n))

    # prepare path list to output
    pretty_paths = []
    for path in finished_paths:
        ppath = Path()
        for net in path:
            if stree.is_sentinel(net):
                net = stree.get_node_from_sentinel(net)
            if stree.is_node(net) or stree.is_trivial_chain(net):
                # if it's a node, add it to the path
                ppath.addNodeHandle(net, stree)
            elif stree.is_chain(net):
                # if it's a chain, we need to write someting like ">Nl>*>Nr"
                nodl = stree.get_bound(net, False, True)
                nodl_s = stree.net_handle_as_string(nodl)
                nodr = stree.get_bound(net, True, False)
                nodr_s = stree.net_handle_as_string(nodr)
                ppath.addNodeHandle(nodl, stree)
                ppath.addNode('*', '>')
                ppath.addNodeHandle(nodr, stree)
        # check if path is mostly traversing nodes in reverse orientation
        if ppath.nreversed() > ppath.size() / 2:
            ppath.flip()
        pretty_paths.append(ppath.print())
    # write the paths out
    outf.write('{}\t{}\n'.format(snarl_id,
                                 ','.join(pretty_paths)))
    npaths += len(pretty_paths)

outf.close()

print('{} paths written in {}.'.format(npaths, args.o))
