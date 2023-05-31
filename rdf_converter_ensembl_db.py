import gzip
import sys
import psutil
import json

input_dir = "./"


def triple(s, p, o):
    print(s, p, o, ".")
    return


def quote(string):
    return "\"" + string + "\""


def strand2faldo(s):
    if s == "1":
        return "faldo:ForwardStrandPosition"
    elif s == "-1":
        return "faldo:ReverseStrandPosition"
    else:
        print(f"Error: Invalid argument \"{s}\" for strand2faldo", file=sys.stderr)
        sys.exit(1)


class Bnode:
    def __init__(self):
        self.properties = []

    def add(self, tpl):
        # `tpl` is a tuple of strings (e.g. ("rdf:type", "owl:Class"))
        self.properties.append(tpl)

    def serialize(self, level=1):
        s = "[\n"
        indent = "    " * level
        for tpl in self.properties:
            s += indent + tpl[0] + " " + tpl[1] + " ;\n"
        s += "    " * (level-1) + "]"
        return s


class Ensembl2turtle:
    prefixes = [
        ['rdf:', '<http://www.w3.org/1999/02/22-rdf-syntax-ns#>'],
        ['rdfs:', '<http://www.w3.org/2000/01/rdf-schema#>'],
        ['faldo:', '<http://biohackathon.org/resource/faldo#>'],
        ['obo:', '<http://purl.obolibrary.org/obo/>'],
        ['dc:', '<http://purl.org/dc/elements/1.1/>'],
        ['dcterms:', '<http://purl.org/dc/terms/>'],
        ['owl:', '<http://www.w3.org/2002/07/owl#>'],
        ['ensg:', '<http://rdf.ebi.ac.uk/resource/ensembl/>'],
        ['terms:', '<http://rdf.ebi.ac.uk/terms/ensembl/>'],
        ['ense:', '<http://rdf.ebi.ac.uk/resource/ensembl.exon/>'],
        ['ensp:', '<http://rdf.ebi.ac.uk/resource/ensembl.protein/>'],
        ['enst:', '<http://rdf.ebi.ac.uk/resource/ensembl.transcript/>'],
        ['ensi:', '<http://identifiers.org/ensembl/>'],
        ['taxonomy:', '<http://identifiers.org/taxonomy/>'],
        ['uniprot:', '<http://purl.uniprot.org/uniprot/>'],
        ['sio:', '<http://semanticscience.org/resource/>'],
        ['skos:', '<http://www.w3.org/2004/02/skos/core#>']
    ]

    def __init__(self, input_dbinfo_file):
        self.flg = True
        self.dbinfo = self.read_dbinfo(input_dbinfo_file)
        self.dbs = self.read_dbs()
        self.taxonomy_id = self.get_taxonomy_id()

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
        for id in gene:
            #print(genes[id])
            sbj = "ensg:" + gene[id][6]
            xref_id = gene[id][4]

            triple(sbj, "a", "terms:EnsemblGene")
            triple(sbj, "terms:biotype", "terms:"+gene[id][0])
            xref_id = gene[id][4]
            if xref_id == "\\N":
                label = gene[id][6]  # Substitute ID for label
            else:
                label = xref[xref_id][2]
            triple(sbj, "rdfs:label", quote(label))
            description = gene[id][5]
            if description == "\\N":
                description = ""
            triple(sbj, "dcterms:description", quote(description))
            triple(sbj, "dcterms:identifier", quote(gene[id][6]))
            triple(sbj, "obo:RO_0002162", "taxonomy:"+self.taxonomy_id)
            # synonym
            synonyms = ", ".join([quote(v[0]) for v in external_synonym.get(xref_id, [])])
            if len(synonyms) >= 1:
                triple(sbj, "skos:altLabel", synonyms)

            # location
            loc = Bnode()
            loc_beg = Bnode()
            loc_end = Bnode()
            loc_beg.add(("a", "faldo:ExactPosition"))
            loc_beg.add(("a", strand2faldo(gene[id][3])))
            loc_beg.add(("faldo:position", gene[id][1]))
            loc_end.add(("a", "faldo:ExactPosition"))
            loc_end.add(("a", strand2faldo(gene[id][3])))
            loc_end.add(("faldo:position", gene[id][2]))
            loc.add(("a", "faldo:Region"))
            loc.add(("faldo:begin", loc_beg.serialize(level=2)))
            loc.add(("faldo:end", loc_end.serialize(level=2)))
            triple(sbj, "faldo:location", loc.serialize())
            # i += 1
            # if i >= 10:
            #     break
        return

    def output_prefixes(self):
        for prefix in Ensembl2turtle.prefixes:
            triple("@prefix", prefix[0], prefix[1])
        return

    def output_turtle(self):
        self.output_prefixes()
        self.rdfize_gene()


def main():
    input_dbinfo_file = sys.argv[1]

    converter = Ensembl2turtle(input_dbinfo_file)
    #converter.rdfize_gene()
    #print(converter.dbs['meta'].keys())
    converter.output_turtle()


if __name__ == '__main__':
    main()
