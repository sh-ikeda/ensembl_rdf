import gzip
import sys
import psutil
import json

input_dir = "./"


def triple(s, p, o):
    print(s, p, o, ".")
    return


def quote(s):
    return "\"" + s + "\""

class Ensembl2turtle:
    prefixes = [
        ['rdfs:', '<http://www.w3.org/2000/01/rdf-schema#>'],
        ['faldo:', '<http://biohackathon.org/resource/faldo#>'],
        ['obo:', '<http://purl.obolibrary.org/obo/>'],
        ['dc:', '<http://purl.org/dc/elements/1.1/>'],
        ['owl:', '<http://www.w3.org/2002/07/owl#>'],
        ['ens_resource:', '<http://rdf.ebi.ac.uk/resource/ensembl/>'],
        ['ens_terms:', '<http://rdf.ebi.ac.uk/terms/ensembl/>'],
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
        taxonomy_ids = [v[2] for k, v in self.dbs["meta"].items() if v[1] == 'species.taxonomy_id']
        if len(taxonomy_ids) >= 2:
            print("Error: `meta` table has multiple taxonomy_id. This seems to be multi-species database.", file=sys.stderr)
            print("taxonomy_ids: ", taxonomy_ids, file=sys.stderr)
            sys.exit(1)
        self.taxonomy_id = taxonomy_ids[0]

    def read_dbinfo(self, input_dbinfo_file):
        dbinfo_dict = {}
        with open(input_dbinfo_file, 'r') as input_dbinfo:
            dbinfo_dict = json.load(input_dbinfo)
            # line = input_dbinfo.readline()
            # while (line):
            #     line = line.rstrip('\n')
            #     sep_line = line.split('\t')

            #     table_name = sep_line[0].replace('.txt.gz', '')
            #     dbinfo_dict[table_name] = {
            #         "key_indices": sep_line[1].split(','),
            #         "val_indices": list(map(int, sep_line[2].split(','))),
            #         "val_names": sep_line[3].split(','),
            #         "table_filename": sep_line[0]
            #     }

            #     line = input_dbinfo.readline()
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
                        print(db, key, vals)
                        self.flg = False
                else:
                    dic[key] = vals
                # for index in val_indices:
                #     dic[key].append(sep_line[index])

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
        #gene_attrib = self.dbs["gene_attrib"]
        #name_id = [k for k, v in self.dbs["attrib_type"].items() if v[0] == 'name'][0]
        #synonym_id = [k for k, v in self.dbs["attrib_type"].items() if v[0] == 'synonym'][0]
        xref = self.dbs["xref"]
        external_synonym = self.dbs["external_synonym"]
        #print(name_id, synonym_id)
        i = 0
        for id in gene:
            #print(genes[id])
            sbj = "ensg:" + gene[id][6]
            triple(sbj, "a", "terms:EnsemblGene")
            triple(sbj, "terms:biotype", "terms:"+gene[id][0])
            xref_id = gene[id][4]
            #triple(sbj, "rdfs:label", gene_attrib[(id, name_id)])
            triple(sbj, "rdfs:label", quote(xref[xref_id][2]))
            triple(sbj, "dcterms:description", quote(gene[id][5]))
            triple(sbj, "dcterms:identifier", quote(gene[id][6]))
            triple(sbj, "obo:RO_0002162", "taxonomy:"+self.taxonomy_id)
            # synonym
            for values in external_synonym[xref_id]:
                triple(sbj, "skos:altLabel", quote(values[0]))
            i += 1
            if i >= 10:
                break
        return


def main():
    input_dbinfo_file = sys.argv[1]

    converter = Ensembl2turtle(input_dbinfo_file)
    converter.rdfize_gene()
    #print(converter.dbs['meta'].keys())


if __name__ == '__main__':
    main()
