import gzip
import sys
import psutil
import json
import rdflib
from rdflib import Graph, Literal

input_dir = "./"


def triple(s, p, o):
    print(s, p, o, ".")
    return


def quote(s):
    return "\"" + s + "\""


class Ensembl2turtle:
    prefixes = [
        ['rdf', 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'],
        ['rdfs', 'http://www.w3.org/2000/01/rdf-schema#'],
        ['faldo', 'http://biohackathon.org/resource/faldo#'],
        ['obo', 'http://purl.obolibrary.org/obo/'],
        ['dc', 'http://purl.org/dc/elements/1.1/'],
        ['dcterms', 'http://purl.org/dc/terms/'],
        ['owl', 'http://www.w3.org/2002/07/owl#'],
        ['ensg', 'http://rdf.ebi.ac.uk/resource/ensembl/'],
        ['terms', 'http://rdf.ebi.ac.uk/terms/ensembl/'],
        ['ense', 'http://rdf.ebi.ac.uk/resource/ensembl.exon/'],
        ['ensp', 'http://rdf.ebi.ac.uk/resource/ensembl.protein/'],
        ['enst', 'http://rdf.ebi.ac.uk/resource/ensembl.transcript/'],
        ['ensi', 'http://identifiers.org/ensembl/'],
        ['taxonomy', 'http://identifiers.org/taxonomy/'],
        ['uniprot', 'http://purl.uniprot.org/uniprot/'],
        ['sio', 'http://semanticscience.org/resource/'],
        ['skos', 'http://www.w3.org/2004/02/skos/core#']
    ]

    def __init__(self, input_dbinfo_file):
        self.flg = True
        self.dbinfo = self.read_dbinfo(input_dbinfo_file)
        self.dbs = self.read_dbs()
        self.taxonomy_id = self.get_taxonomy_id()
        self.namespace_dict = {}
        self.initialize_namespace_dict()
        self.graph = Graph()
        self.initialize_graph()

    # self.namespace_dict = {
    #  "ensg": Namespace("http://rdf.ebi.ac.uk/resource/ensembl/"), ..}
    def initialize_namespace_dict(self):
        for prefix in Ensembl2turtle.prefixes:
            ns = rdflib.Namespace(prefix[1])
            self.namespace_dict[prefix[0]] = ns
        return

    def initialize_graph(self):
        for abbrev, ns in self.namespace_dict.items():
            self.graph.namespace_manager.bind(abbrev, ns)
        return

    def get_taxonomy_id(self):
        taxonomy_ids = [v[2] for k, v in self.dbs["meta"].items() if v[1] == 'species.taxonomy_id']
        if len(taxonomy_ids) >= 2:
            print("Error: `meta` table has multiple taxonomy_id. This seems to be multi-species database.", file=sys.stderr)
            print("taxonomy_ids: ", taxonomy_ids, file=sys.stderr)
            sys.exit(1)

        return taxonomy_ids[0]

    def read_dbinfo(self, input_dbinfo_file):
        dbinfo_dict = {}
        with open(input_dbinfo_file, 'r') as input_dbinfo:
            dbinfo_dict = json.load(input_dbinfo)
        return dbinfo_dict

    def read_db(self, db):
        print("Reading DB: ", db, file=sys.stderr)
        dic = {}
        table_file = self.dbinfo[db]["filename"]
        key_indices = self.dbinfo[db]["key_indices"]
        val_indices = [v["index"] for v in self.dbinfo[db]["values"]]
        with gzip.open(table_file, 'rt') as input_table:
            line = input_table.readline()
            while (line):
                line = line.rstrip('\n')
                sep_line = line.split('\t')
                # gene_attrib の sep_line: [gene_id, attrib_type_id, value]
                key_list = [sep_line[i] for i in key_indices]
                if len(key_list) >= 2:
                    key = tuple(key_list)
                else:
                    key = key_list[0]
                vals = [sep_line[i] for i in val_indices]
                if self.dbinfo[db].get("list", False):
                    if key in dic:
                        dic[key].append(vals)
                    else:
                        dic[key] = [vals]
                    if self.flg:
                        print(db, key, vals, file=sys.stderr)
                        self.flg = False
                else:
                    dic[key] = vals

                line = input_table.readline()

        return dic

    def read_dbs(self):
        db_dics = {}
        for db in self.dbinfo:
            db_dics[db] = self.read_db(db)
            process = psutil.Process()
            memory_usage = process.memory_info().rss  # バイト単位でのメモリ使用量
            # print(db)
            # print(f'memory: {memory_usage / (1024*1024)} MB')

        return db_dics

    def rdfize_gene(self):
        gene = self.dbs["gene"]
        xref = self.dbs["xref"]
        external_synonym = self.dbs["external_synonym"]
        #print(name_id, synonym_id)
        i = 0
        nsd = self.namespace_dict
        for id in gene:
            #print(genes[id])
            sbj = nsd["ensg"][gene[id][6]]
            xref_id = gene[id][4]

            self.graph.add((sbj, nsd["rdf"]["type"], nsd["terms"]["EnsemblGene"]))
            self.graph.add((sbj, nsd["terms"]["biotype"], nsd["terms"][gene[id][0]]))
            self.graph.add((sbj, nsd["rdfs"]["label"], Literal(xref[xref_id][2])))
            self.graph.add((sbj, nsd["rdfs"]["label"], Literal(xref[xref_id][2])))
            self.graph.add((sbj, nsd["dcterms"]["description"], Literal(gene[id][5])))
            self.graph.add((sbj, nsd["dcterms"]["identifier"], Literal(gene[id][6])))
            self.graph.add((sbj, nsd["obo"]["RO_0002162"], nsd["taxonomy"][self.taxonomy_id]))
            # synonym
            for syn in external_synonym[xref_id]:
                self.graph.add((sbj, nsd["skos"]["altLabel"], Literal(syn[0])))

            # sbj = "ensg:" + gene[id][6]
            # triple(sbj, "a", "terms:EnsemblGene")
            # triple(sbj, "terms:biotype", "terms:"+gene[id][0])
            # xref_id = gene[id][4]
            # triple(sbj, "rdfs:label", quote(xref[xref_id][2]))
            # triple(sbj, "dcterms:description", quote(gene[id][5]))
            # triple(sbj, "dcterms:identifier", quote(gene[id][6]))
            # triple(sbj, "obo:RO_0002162", "taxonomy:"+self.taxonomy_id)
            # # synonym
            # synonyms = ", ".join([quote(v[0]) for v in external_synonym[xref_id]])
            # triple(sbj, "skos:altLabel", synonyms)
            i += 1
            if i >= 10:
                break
        return

    def output_graph(self):
        print(self.graph.serialize(format="turtle"))
        return


def main():
    input_dbinfo_file = sys.argv[1]

    converter = Ensembl2turtle(input_dbinfo_file)
    converter.rdfize_gene()
    #print(converter.dbs['meta'].keys())
    converter.output_graph()


if __name__ == '__main__':
    main()
