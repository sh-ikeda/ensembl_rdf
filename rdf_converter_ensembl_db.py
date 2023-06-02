import gzip
import sys
import psutil
import json
import re
import datetime

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
        ['so:', '<http://purl.obolibrary.org/obo/so#>'],
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

    # Keys are `code` of the attrib_type table to be used as transcript flags
    transcript_flags = {
        "gencode_basic": {"is_bool": True, "term": "terms:GENCODEBasic"},
        "appris": {"is_bool": False, "term": "terms:APPRIS"},
        "TSL": {"is_bool": False, "term": "terms:TSL"},
        "is_canonical": {"is_bool": True, "term": "terms:EnsemblCanonical"},
        "MANE_Select": {"is_bool": True, "term": "terms:MANESelect"},
        "mRNA_start_NF": {"is_bool": False, "term": "terms:Incomplete"},
        "mRNA_end_NF": {"is_bool": False, "term": "terms:Incomplete"},
        "cds_start_NF": {"is_bool": False, "term": "terms:Incomplete"},
        "cds_end_NF": {"is_bool": False, "term": "terms:Incomplete"}
    }

    def __init__(self, input_dbinfo_file):
        self.flg = True
        self.debug = False
        self.dbinfo = self.load_dbinfo(input_dbinfo_file)
        self.dbs = self.load_dbs()
        self.taxonomy_id = self.get_taxonomy_id()
        self.ensembl_version = self.get_ensembl_version()
        self.production_name = self.get_production_name()
        self.flag_dic = {}
        self.init_flags()

    # self.flag_dic は transcript_flags のキーを attrib_type_id で置換したもの
    # ただし、transcript_flags のキーは "name" として保存
    def init_flags(self):
        attrib_type = self.dbs["attrib_type"]
        # {"428": ["TSL", "Transcript Support Level", "(desc)"]}
        for id, vals in attrib_type.items():
            if vals[0] in Ensembl2turtle.transcript_flags:
                self.flag_dic[id] = Ensembl2turtle.transcript_flags[vals[0]]
                self.flag_dic[id]["name"] = vals[0]
        return

    def get_ensembl_version(self):
        ensembl_version = [v[2] for k, v in self.dbs["meta"].items() if v[1] == 'schema_version']
        return ensembl_version[0]

    def get_taxonomy_id(self):
        taxonomy_ids = [v[2] for k, v in self.dbs["meta"].items() if v[1] == 'species.taxonomy_id']
        if len(taxonomy_ids) >= 2:
            print("Error: `meta` table has multiple taxonomy_id. This seems to be multi-species database.", file=sys.stderr)
            print("taxonomy_ids: ", taxonomy_ids, file=sys.stderr)
            sys.exit(1)

        return taxonomy_ids[0]

    def get_production_name(self):
        production_names = [v[2] for k, v in self.dbs["meta"].items() if v[1] == 'species.production_name']
        return production_names[0]

    def load_dbinfo(self, input_dbinfo_file):
        dbinfo_dict = {}
        with open(input_dbinfo_file, 'r') as input_dbinfo:
            dbinfo_dict = json.load(input_dbinfo)
        return dbinfo_dict

    def load_db(self, db):
        dt_now = datetime.datetime.now()
        print(f"[{dt_now}] Loading DB: {db}", file=sys.stderr)
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
                    # if self.flg:
                    #     print(db, key, vals, file=sys.stderr)
                    #     self.flg = False
                else:
                    dic[key] = vals

                line = input_table.readline()

        return dic

    def load_dbs(self):
        db_dics = {}
        for db in self.dbinfo:
            db_dics[db] = self.load_db(db)
            process = psutil.Process()
            memory_usage = process.memory_info().rss  # バイト単位でのメモリ使用量
            print(f'memory: {memory_usage / (1024*1024)} MB', file=sys.stderr)

        return db_dics

    def rdfize_gene(self):
        gene = self.dbs["gene"]
        xref = self.dbs["xref"]
        seq_region = self.dbs["seq_region"]
        external_synonym = self.dbs["external_synonym"]
        i = 0
        for id in gene:
            sbj = "ensg:" + gene[id][6]
            xref_id = gene[id][4]

            triple(sbj, "a", "terms:EnsemblGene")
            triple(sbj, "terms:biotype", "terms:"+gene[id][0])
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
            location = self.create_location_str(gene[id][1],
                                                gene[id][2],
                                                gene[id][3],
                                                gene[id][7])
            triple(sbj, "faldo:location", location)
            i += 1
            if self.debug and i >= 10:
                break
        return

    def rdfize_transcript(self):
        transcript = self.dbs["transcript"]
        transcript_attrib = self.dbs["transcript_attrib"]
        xref = self.dbs["xref"]
        gene = self.dbs["gene"]
        translation = self.dbs["translation"]
        i = 0
        for id in transcript:
            sbj = "enst:" + transcript[id][7]
            xref_id = transcript[id][4]

            triple(sbj, "a", "terms:EnsemblTranscript")
            triple(sbj, "terms:biotype", "terms:"+transcript[id][5])
            if xref_id == "\\N":
                label = transcript[id][7]  # Substitute ID for label
            else:
                label = xref[xref_id][2]
            triple(sbj, "rdfs:label", quote(label))
            triple(sbj, "dcterms:identifier", quote(transcript[id][7]))
            triple(sbj, "so:transcribed_from", "ensg:"+gene[transcript[id][0]][6])
            translates_to = transcript[id][6]
            if translates_to != "\\N":
                triple(sbj, "so:translates_to", "ensp:"+translation[translates_to][1])

            # location
            location = self.create_location_str(transcript[id][1],
                                                transcript[id][2],
                                                transcript[id][3],
                                                transcript[id][8])
            triple(sbj, "faldo:location", location)

            # flag
            attribs = transcript_attrib.get(id, [])
            for attrib in attribs:
                flag = Bnode()
                flag.add(("a", "terms:TranscriptFlag"))
                if attrib[0] in self.flag_dic:
                    flag.add(("a", self.flag_dic[attrib[0]]["term"]))
                    if self.flag_dic[attrib[0]]["is_bool"]:
                        flag.add(("sio:SIO_000300", quote("true")))
                    elif self.flag_dic[attrib[0]]["name"] == "TSL":
                        tsl = re.sub(r" .*", "", attrib[1])
                        comment = ""
                        match = re.search(r'\((.*?)\)', attrib[1])
                        if match:
                            comment = match.group(1)
                        flag.add(("sio:SIO_000300", quote(tsl)))
                        if len(comment) > 0:
                            flag.add(("rdfs:comment", quote(comment)))
                    elif self.flag_dic[attrib[0]]["term"] == "terms:Incomplete":
                        flag.add(("sio:SIO_000300", quote(self.flag_dic[attrib[0]]["name"])))
                    else:
                        flag.add(("sio:SIO_000300", quote(attrib[1])))
                    triple(sbj, "terms:has_transcript_flag", flag.serialize())

            i += 1
            if self.debug and i >= 10:
                break
        return

    def create_location_str(self, beg, end, strand, seq_region_id, level=1):
        loc = Bnode()
        loc_beg = Bnode()
        loc_end = Bnode()

        seq_region = self.dbs["seq_region"]
        coord_system = self.dbs["coord_system"]
        chromosome_name = seq_region[seq_region_id][0]
        coord_system_id = seq_region[seq_region_id][1]
        # e.g. "GRCm38"
        coord_system_version = coord_system[coord_system_id][2]
        # e.g. <http://rdf.ebi.ac.uk/resource/ensembl/109/mus_musculus/GRCm38/Y>
        chromosome_url = "<http://rdf.ebi.ac.uk/resource/ensembl/"+self.ensembl_version+"/"+self.production_name+"/"+coord_system_version+"/"+chromosome_name+">"
        # For LRG, <http://rdf.ebi.ac.uk/resource/ensembl/109/homo_sapiens/LRG_1>">"
        if coord_system[coord_system_id][1] == "lrg":
            chromosome_url = "<http://rdf.ebi.ac.uk/resource/ensembl/"+self.ensembl_version+"/"+self.production_name+"/"+chromosome_name+">"

        loc_beg.add(("a", "faldo:ExactPosition"))
        loc_beg.add(("a", strand2faldo(strand)))
        loc_beg.add(("faldo:position", beg))
        loc_beg.add(("faldo:reference", chromosome_url))
        loc_end.add(("a", "faldo:ExactPosition"))
        loc_end.add(("a", strand2faldo(strand)))
        loc_end.add(("faldo:position", end))
        loc_end.add(("faldo:reference", chromosome_url))
        loc.add(("a", "faldo:Region"))
        loc.add(("faldo:begin", loc_beg.serialize(level=level+1)))
        loc.add(("faldo:end", loc_end.serialize(level=level+1)))
        return loc.serialize()

    def rdfize_translation(self):
        transcript = self.dbs["transcript"]
        xref = self.dbs["xref"]
        translation = self.dbs["translation"]
        i = 0
        for id in translation:
            sbj = "ensp:" + translation[id][1]

            triple(sbj, "a", "terms:EnsemblProtein")
            triple(sbj, "dcterms:identifier", quote(translation[id][1]))
            triple(sbj, "so:translation_of", "enst:"+transcript[translation[id][0]][7])
            i += 1
            if self.debug and i >= 10:
                break
        return

    def rdfize_exon(self):
        transcript = self.dbs["transcript"]
        exon = self.dbs["exon"]
        translation = self.dbs["translation"]
        i = 0
        for id in exon:
            sbj = "ense:" + exon[id][3]

            triple(sbj, "a", "terms:EnsemblExon")
            triple(sbj, "a", "obo:SO_0000147")
            triple(sbj, "dcterms:identifier", quote(exon[id][3]))

            # location
            location = self.create_location_str(exon[id][0],
                                                exon[id][1],
                                                exon[id][2],
                                                exon[id][4])
            triple(sbj, "faldo:location", location)
            i += 1
            if self.debug and i >= 10:
                break
        return

    def rdfize_exon_transcript(self):
        exon_transcript = self.dbs["exon_transcript"]
        transcript = self.dbs["transcript"]
        exon = self.dbs["exon"]
        i = 0
        for id in exon_transcript:
            exon_id = id[0]
            transcript_id = id[1]
            exon_stable_id = exon[exon_id][3]
            transcript_stable_id = transcript[transcript_id][7]
            rank = exon_transcript[id][0]
            ordered_exon_uri = "<http://rdf.ebi.ac.uk/resource/ensembl.transcript/"+transcript_stable_id+"#Exon_"+rank+">"
            exon_uri = "ense:" + exon_stable_id
            transcript_uri = "enst:" + transcript_stable_id

            triple(ordered_exon_uri, "a", "terms:EnsemblOrderedExon")
            triple(ordered_exon_uri, "a", "sio:SIO_001261")
            triple(ordered_exon_uri, "sio:SIO_000628", exon_uri)
            triple(ordered_exon_uri, "sio:SIO_000300", rank)

            triple(transcript_uri, "so:has_part", exon_uri)
            triple(transcript_uri, "sio:SIO_000974", ordered_exon_uri)
            i += 1
            if self.debug and i >= 10:
                break
        return

    def rdfize_xref(self):
        gene = self.dbs["gene"]
        transcript = self.dbs["transcript"]
        translation = self.dbs["translation"]
        xref = self.dbs["xref"]
        object_xref = self.dbs["object_xref"]
        external_db = self.dbs["external_db"]
        i = 0
        for id in object_xref:
            xref_id = object_xref[id][2]
            target_id = object_xref[id][0]
            if object_xref[id][1] == "Gene":
                target_url = "ensg:" + gene[target_id][6]
            elif object_xref[id][1] == "Transcript":
                target_url = "enst:" + transcript[target_id][7]
            elif object_xref[id][1] == "Translation":
                target_url = "ensp:" + translation[target_id][1]
            else:
                continue
            xref_node = Bnode()
            xref_node.add(("terms:id_of", quote(external_db[xref[xref_id][0]][0])))  # FIXME
            xref_node.add(("dcterms:identifier", quote(xref[xref_id][1])))
            triple(target_url, "rdfs:seeAlso", xref_node.serialize())
            i += 1
            if self.debug and i >= 10:
                break
        return

    def output_prefixes(self):
        for prefix in Ensembl2turtle.prefixes:
            triple("@prefix", prefix[0], prefix[1])
        return

    def output_turtle(self):
        self.output_prefixes()
        dt_now = datetime.datetime.now()
        print(f"[{dt_now}] Output turtle: gene", file=sys.stderr)
        self.rdfize_gene()
        dt_now = datetime.datetime.now()
        print(f"[{dt_now}] Output turtle: transcript", file=sys.stderr)
        self.rdfize_transcript()
        dt_now = datetime.datetime.now()
        print(f"[{dt_now}] Output turtle: translation", file=sys.stderr)
        self.rdfize_translation()
        dt_now = datetime.datetime.now()
        print(f"[{dt_now}] Output turtle: exon", file=sys.stderr)
        self.rdfize_exon()
        dt_now = datetime.datetime.now()
        print(f"[{dt_now}] Output turtle: exon_transcript", file=sys.stderr)
        self.rdfize_exon_transcript()
        dt_now = datetime.datetime.now()
        print(f"[{dt_now}] Output turtle: xref", file=sys.stderr)
        self.rdfize_xref()
        dt_now = datetime.datetime.now()
        print(f"[{dt_now}] Done.", file=sys.stderr)


def main():
    input_dbinfo_file = sys.argv[1]

    converter = Ensembl2turtle(input_dbinfo_file)
    #converter.rdfize_gene()
    #print(converter.dbs['meta'].keys())
    converter.output_turtle()


if __name__ == '__main__':
    main()
