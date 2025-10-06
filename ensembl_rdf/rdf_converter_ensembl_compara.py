import os
import sys
from utils import log_time, percent_encode, quote_str
from ensembl2turtle import Ensembl2turtle


class Compara2turtle(Ensembl2turtle):
    def __init__(self, input_dbinfo_file, input_data_dir):
        super().__init__(input_dbinfo_file, input_data_dir)
        self.prefixes += [
            ['genetree:', '<http://identifiers.org/genetree/>'],
            ['genetree_plants:', '<http://identifiers.org/genetree.plants/>'],
            ['genetree_pan:', '<http://identifiers.org/genetree.bacteria/>'],
            ['genetree_protists:', '<http://identifiers.org/genetree.protists/>'],
            ['genetree_metazoa:', '<http://identifiers.org/genetree.metazoa/>'],
            ['genetree_fungi:', '<http://identifiers.org/genetree.fungi/>'],
            ['orth:', '<http://purl.org/net/orth#>']
        ]

    def rdfize_genetree(self):
        gene_tree_node = self.dbs["gene_tree_node"]
        gene_tree_root = self.dbs["gene_tree_root"]
        seq_member = self.dbs["seq_member"]
        f = open(os.path.join(self.input_data_dir, "genetree.ttl"), mode="w")
        self.output_file = f
        self.output_prefixes()
        for id in gene_tree_root:
            tree_type = gene_tree_root[id][1]
            gene_tree_id = gene_tree_root[id][2]
            if gene_tree_id == "\\N" or tree_type != "tree":
                continue
            root_prefix = "genetree:"
            root_uri = root_prefix + gene_tree_id
            self.triple(root_uri, "a", "orth:OrthologsCluster")
            self.triple(root_uri, "dcterms:identifier", quote_str(gene_tree_id))

        for id in gene_tree_node:
            root_id = gene_tree_node[id][0]
            gene_tree_id = gene_tree_root[root_id][2]
            seq_member_id = gene_tree_node[id][1]
            if seq_member_id == "\\N":  # not a sequence node
                continue
            ensembl_id = seq_member[seq_member_id][0]
            member_type = gene_tree_root[root_id][0]
            if member_type == "protein":
                node_prefix = "ensp:"
            elif member_type == "ncrna":
                node_prefix = "enst:"
            node_uri = node_prefix + ensembl_id
            root_prefix = "genetree:"
            root_uri = root_prefix + gene_tree_id
            self.triple(root_uri, "orth:hasHomologousMember", node_uri)

        return

    def output_turtle(self):
        log_time(f"Output turtle: genetree")
        self.rdfize_genetree()
        return

def main():
    input_dbinfo_file = sys.argv[1]
    input_data_dir = sys.argv[2]

    converter = Compara2turtle(input_dbinfo_file, input_data_dir)
    converter.output_turtle()


if __name__ == '__main__':
    main()
