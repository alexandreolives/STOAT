import argparse
import os
import time
import bdsg
from bdsg.bdsg import PackedGraph

class Graph_chunk :
    def __init__(self, path_pangenome, path_distance) -> None:
        self.list_idx_decompose = []
        self.pg = self._initilise_pangenome(path_pangenome)
        self.stree = self._initilise_distance(path_distance)
        self.root = self.stree.get_root()

    def _initilise_pangenome(self, path_pangenome) :
        pg = bdsg.bdsg.PackedGraph()
        pg.deserialize(path_pangenome)
        return pg

    def _initilise_distance(self, path_distance) :
        stree = bdsg.bdsg.SnarlDistanceIndex()
        stree.deserialize(path_distance)
        return stree
    
    def save_snarl_tree_node(self, net, snarls) -> bool :
        if self.stree.is_snarl(net):
            snarls.append(net)
        if not self.stree.is_node(net) and self.stree.is_sentinel(net):
            self.stree.for_each_child(net, self.save_snarl_tree_node)
        return (True)
    
    def loop_over_children(self) -> list :
        snarls = []
        self.stree.for_each_child(self.root, self.save_snarl_tree_node)
        snarls_length = len(snarls)
        print('{} snarls found'.format(snarls_length))

        return snarls

    def indentify_chunk(self, snarls, distance_max=100000) -> list :
        curr_dist = 0
        # for each snarl, lists paths through the netgraph
        for _, snarl in enumerate(snarls):

            # create a snarl ID as LEFT_RIGTH bondary nodes
            sstart = self.stree.get_bound(snarl, False, True)
            sstart = self.stree.get_node_from_sentinel(sstart)
            send = self.stree.get_bound(snarl, True, True)
            send = self.stree.get_node_from_sentinel(send)
            curr_dist += send - sstart

            if curr_dist > distance_max :
                curr_dist = 0
                self.list_idx_decompose.append(save_snarl)
            
            save_snarl = send

    def decompose(self, directory='output_pangenome_decomposed') :
        # Ensure the directory exists
        if not os.path.exists(directory):
            os.makedirs(directory)

        for index, sublist in enumerate(self.list_idx_decompose):
            file_path = os.path.join(directory, f"sublist_{index + 1}.txt")
            
            with open(file_path, 'w') as file:
                for item in sublist:
                    file.write(f"{item}\n")

if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description="Parse and analyse snarl from vcf file")
    parser.add_argument('-p', "--pangenome", help='the input pangenome .pg file', required=True)
    parser.add_argument('-d', "--distance", help='the input distance index .dist file', required=True)
    parser.add_argument('-o', "--output", help='the output dir where all .dist and .pg will be decompose', required=False)

    args = parser.parse_args()
    graph = Graph_chunk()

    snarls = graph.loop_over_children()
    graph.indentify_chunk()

    if args.o :
        graph.decompose(args.o)
    else :
        graph.decompose()


