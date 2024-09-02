import argparse
import os
import time
import bdsg
from bdsg.bdsg import PackedGraph

class Graph_chunk :
    def __init__(self, pg_file_path,  dist_file_path) -> None:

        self.pg = PackedGraph()
        self.pg.deserialize(pg_file_path)

        self.dist = PackedGraph()
        self.dist.deserialize(dist_file_path) 

    def decompose(self) :
        ...


if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description="Parse and analyse snarl from vcf file")
    parser.add_argument('-p', "--pangenome", help='the input pangenome .pg file', required=True)
    parser.add_argument('-d', "--distance", help='the input distance index .dist file',required=True)
    args = parser.parse_args()
    #graph = Graph_chunk(args.pangenome)

    pg = bdsg.bdsg.PackedGraph()
    pg.deserialize(args.pangenome)

    stree = bdsg.bdsg.SnarlDistanceIndex()
    stree.deserialize(args.distance)
    root = stree.get_root()

    snarls = []
    sentinels = []
    nodes = []
    itr = 0
    def save_snarl_tree_node(net):
        if len(snarls) >= 20 :
            return
        if stree.is_sentinel(net):
            sentinels.append(net)

        if stree.is_snarl(net):
            snarls.append(net)

        if stree.is_node(net):
            nodes.append(net)

        if not stree.is_node(net) and not stree.is_sentinel(net):
            stree.for_each_child(net, save_snarl_tree_node)
        return (True)

    stree.for_each_child(root, save_snarl_tree_node)
    snarls_length = len(snarls)
    print("snarl : ", snarls)
    print("sentinels : ", sentinels)
    print("nodes : ", nodes)

    # for each snarl, lists paths through the netgraph and write to output TSV
    for idx, snarl in enumerate(snarls):

        # create a snarl ID as LEFT_RIGTH bondary nodes
        sstart = stree.get_bound(snarl, False, True)
        sstart = stree.get_node_from_sentinel(sstart)
        send = stree.get_bound(snarl, True, True)
        send = stree.get_node_from_sentinel(send)
        snarl_id = '{}_{}'.format(stree.node_id(sstart),
                                stree.node_id(send))
        print("snarl_id : ", snarl_id)

