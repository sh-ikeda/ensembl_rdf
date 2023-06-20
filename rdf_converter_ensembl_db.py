import gzip
import sys
import os
import psutil
import json
import re
import datetime

input_dir = "./"


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
        "gencode_basic": {"GENCODE basic": "terms:GENCODEBasic"},
        "appris": {
            "principal1": "terms:Principal1",
            "principal2": "terms:Principal2",
            "principal3": "terms:Principal3",
            "principal4": "terms:Principal4",
            "principal5": "terms:Principal5",
            "alternative1": "terms:Alternative1",
            "alternative2": "terms:Alternative2"
        },
        "TSL": {
            "tsl1": "terms:TSL1",
            "tsl2": "terms:TSL2",
            "tsl3": "terms:TSL3",
            "tsl4": "terms:TSL4",
            "tsl5": "terms:TSL5",
            "tslNA": "terms:TSLNA"
        },
        "is_canonical": {"1": "terms:EnsemblCanonical"},
        "remark": {"MANE_select": "terms:MANESelect"},
        "mRNA_start_NF": {"1": "terms:Incomplete"},
        "mRNA_end_NF": {"1": "terms:Incomplete"},
        "cds_start_NF": {"1": "terms:Incomplete"},
        "cds_end_NF": {"1": "terms:Incomplete"}
    }

    hco_chr_names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                     "11", "12", "13", "14", "15", "16", "17", "18",
                     "19", "20", "21", "22", "X", "Y", "MT"]

    def __init__(self, input_dbinfo_file):
        self.dbinfo = self.load_dbinfo(input_dbinfo_file)
        self.dbs = self.load_dbs()
        self.taxonomy_id = self.get_taxonomy_id()
        self.ensembl_version = self.get_ensembl_version()
        self.production_name = self.get_production_name()
        self.xref_url_dic = {}
        self.init_xref_url_dic()
        self.xrefed_dbs = {"Gene": {}, "Transcript": {}, "Translation": {}}
        self.not_xrefed_dbs = {"Gene": {}, "Transcript": {}, "Translation": {}}
        self.output_file = sys.stdout

    def init_xref_url_dic(self):
        BASE_DIR = os.path.dirname(os.path.abspath(__file__)) + "/"
        XREF_URL_DIC_TSV = "external_db_url.tsv"

        with open(BASE_DIR+XREF_URL_DIC_TSV, "r") as input_table:
            line = input_table.readline()
            while (line):
                line = line.rstrip('\n')
                sep_line = line.split('\t')
                if sep_line[1] != "":
                    self.xref_url_dic[sep_line[0]] = sep_line[1]
                line = input_table.readline()
        return

    def triple(self, s, p, o):
        print(s, p, o, ".", file=self.output_file)
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
        f = open("gene.ttl", mode="w")
        self.output_file = f
        self.output_prefixes()
        for id in gene:
            sbj = "ensg:" + gene[id][6]
            xref_id = gene[id][4]

            self.triple(sbj, "a", "terms:EnsemblGene")
            self.triple(sbj, "a", "terms:"+gene[id][0])
            self.triple(sbj, "terms:biotype", "terms:"+gene[id][0])
            if xref_id == "\\N":
                label = gene[id][6]  # Substitute ID for label
            else:
                label = xref[xref_id][2]
            self.triple(sbj, "rdfs:label", quote(label))
            description = gene[id][5]
            if description == "\\N":
                description = ""
            self.triple(sbj, "dcterms:description", quote(description))
            self.triple(sbj, "dcterms:identifier", quote(gene[id][6]))
            self.triple(sbj, "obo:RO_0002162", "taxonomy:"+self.taxonomy_id)
            # synonym
            synonyms = ", ".join([quote(v[0]) for v in external_synonym.get(xref_id, [])])
            if len(synonyms) >= 1:
                self.triple(sbj, "skos:altLabel", synonyms)

            # location
            seq_region_id = gene[id][7]
            chromosome_urls = self.seq_region_id_to_chr(seq_region_id)
            location = self.create_location_str(gene[id][1],
                                                gene[id][2],
                                                gene[id][3],
                                                chromosome_urls)
            self.triple(sbj, "faldo:location", location)
            for chromosome_url in chromosome_urls:
                self.triple(sbj, "so:part_of", chromosome_url)
        self.output_file = sys.stdout
        f.close()
        return

    def rdfize_transcript(self):
        transcript = self.dbs["transcript"]
        transcript_attrib = self.dbs["transcript_attrib"]
        xref = self.dbs["xref"]
        gene = self.dbs["gene"]
        translation = self.dbs["translation"]
        attrib_type = self.dbs["attrib_type"]
        f = open("transcript.ttl", mode="w")
        self.output_file = f
        self.output_prefixes()
        for id in transcript:
            sbj = "enst:" + transcript[id][7]
            xref_id = transcript[id][4]

            self.triple(sbj, "a", "terms:EnsemblTranscript")
            self.triple(sbj, "a", "terms:"+transcript[id][5])
            self.triple(sbj, "terms:biotype", "terms:"+transcript[id][5])
            if xref_id == "\\N":
                label = transcript[id][7]  # Substitute ID for label
            else:
                label = xref[xref_id][2]
            self.triple(sbj, "rdfs:label", quote(label))
            self.triple(sbj, "dcterms:identifier", quote(transcript[id][7]))
            self.triple(sbj, "so:transcribed_from", "ensg:"+gene[transcript[id][0]][6])
            translates_to = transcript[id][6]
            if translates_to != "\\N":
                self.triple(sbj, "so:translates_to", "ensp:"+translation[translates_to][1])

            # location
            chromosome_urls = self.seq_region_id_to_chr(transcript[id][8])
            location = self.create_location_str(transcript[id][1],
                                                transcript[id][2],
                                                transcript[id][3],
                                                chromosome_urls)
            self.triple(sbj, "faldo:location", location)

            # flag
            attribs = transcript_attrib.get(id, [])
            flag_dic = Ensembl2turtle.transcript_flags
            for attrib in attribs:
                attrib_code = attrib_type[attrib[0]][0]
                if attrib_code in flag_dic:
                    attrib_val = attrib[1]
                    if attrib_code == "TSL":
                        comment = ""
                        match = re.search(r'\((.*?)\)', attrib_val)
                        attrib_val = re.sub(r" .*", "", attrib[1])
                        if match:
                            comment = match.group(1)
                            statement = "<http://rdf.ebi.ac.uk/resource/ensembl.transcript/#_" + transcript[id][7] + "-has_transcript_flag-"+attrib_val+">"
                            self.triple(statement, "a", "rdf:Statement")
                            self.triple(statement, "rdf:subject", sbj)
                            self.triple(statement, "rdf:predicate", "terms:has_transcript_flag")
                            self.triple(statement, "rdf:object", flag_dic[attrib_code][attrib_val])
                            self.triple(statement, "rdfs:comment", quote(comment))
                    elif attrib_code == "remark":
                        if attrib_val != "MANE_select":
                            continue
                    self.triple(sbj, "terms:has_transcript_flag", flag_dic[attrib_code][attrib_val])
                    # try:
                    #     self.triple(sbj, "terms:has_transcript_flag", flag_dic[attrib_code][attrib_val])
                    # except KeyError as e:
                    #     print(sbj, attrib_code, attrib_val, file=sys.stderr)
                    #     sys.exit()
        self.output_file = sys.stdout
        f.close()
        return

    def seq_region_id_to_chr(self, seq_region_id):
        seq_region = self.dbs["seq_region"]
        coord_system = self.dbs["coord_system"]
        chromosome_name = seq_region[seq_region_id][0]
        coord_system_id = seq_region[seq_region_id][1]
        chromosome_urls = []
        # e.g. "GRCm38"
        coord_system_version = coord_system[coord_system_id][2]
        # e.g. <http://rdf.ebi.ac.uk/resource/ensembl/109/mus_musculus/GRCm38/Y>
        chromosome_url = "<http://rdf.ebi.ac.uk/resource/ensembl/"+self.ensembl_version+"/"+self.production_name+"/"+coord_system_version+"/"+chromosome_name+">"
        # For LRG, <http://rdf.ebi.ac.uk/resource/ensembl/109/homo_sapiens/LRG_1>">"
        if coord_system[coord_system_id][1] == "lrg":
            chromosome_url = "<http://rdf.ebi.ac.uk/resource/ensembl/"+self.ensembl_version+"/"+self.production_name+"/"+chromosome_name+">"
        chromosome_urls.append(chromosome_url)

        if self.taxonomy_id == "9606":
            if chromosome_name in Ensembl2turtle.hco_chr_names:
                hco_url = "<http://identifiers.org/hco/"+chromosome_name+"/GRCh38>"
                chromosome_urls.append(hco_url)

        return chromosome_urls

    def create_location_str(self, beg, end, strand, chromosome_urls, level=1):
        loc = Bnode()
        loc_beg = Bnode()
        loc_end = Bnode()

        loc_beg.add(("a", "faldo:ExactPosition"))
        loc_beg.add(("a", strand2faldo(strand)))
        loc_beg.add(("faldo:position", beg))
        loc_end.add(("a", "faldo:ExactPosition"))
        loc_end.add(("a", strand2faldo(strand)))
        loc_end.add(("faldo:position", end))

        for chromosome_url in chromosome_urls:
            loc_beg.add(("faldo:reference", chromosome_url))
            loc_end.add(("faldo:reference", chromosome_url))

        loc.add(("a", "faldo:Region"))
        loc.add(("faldo:begin", loc_beg.serialize(level=level+1)))
        loc.add(("faldo:end", loc_end.serialize(level=level+1)))
        return loc.serialize()

    def rdfize_translation(self):
        transcript = self.dbs["transcript"]
        xref = self.dbs["xref"]
        translation = self.dbs["translation"]
        f = open("translation.ttl", mode="w")
        self.output_file = f
        self.output_prefixes()
        for id in translation:
            sbj = "ensp:" + translation[id][1]

            self.triple(sbj, "a", "terms:EnsemblProtein")
            self.triple(sbj, "dcterms:identifier", quote(translation[id][1]))
            self.triple(sbj, "so:translation_of", "enst:"+transcript[translation[id][0]][7])
        self.output_file = sys.stdout
        f.close()
        return

    def rdfize_exon(self):
        transcript = self.dbs["transcript"]
        exon = self.dbs["exon"]
        translation = self.dbs["translation"]
        f = open("exon.ttl", mode="w")
        self.output_file = f
        self.output_prefixes()
        for id in exon:
            sbj = "ense:" + exon[id][3]

            self.triple(sbj, "a", "terms:EnsemblExon")
            self.triple(sbj, "a", "obo:SO_0000147")
            self.triple(sbj, "dcterms:identifier", quote(exon[id][3]))

            # location
            chromosome_urls = self.seq_region_id_to_chr(exon[id][4])
            location = self.create_location_str(exon[id][0],
                                                exon[id][1],
                                                exon[id][2],
                                                chromosome_urls)
            self.triple(sbj, "faldo:location", location)
        self.output_file = sys.stdout
        f.close()
        return

    def rdfize_exon_transcript(self):
        exon_transcript = self.dbs["exon_transcript"]
        transcript = self.dbs["transcript"]
        exon = self.dbs["exon"]
        f = open("exon_transcript.ttl", mode="w")
        self.output_file = f
        self.output_prefixes()
        for id in exon_transcript:
            exon_id = id[0]
            transcript_id = id[1]
            exon_stable_id = exon[exon_id][3]
            transcript_stable_id = transcript[transcript_id][7]
            rank = exon_transcript[id][0]
            ordered_exon_uri = "<http://rdf.ebi.ac.uk/resource/ensembl.transcript/"+transcript_stable_id+"#Exon_"+rank+">"
            exon_uri = "ense:" + exon_stable_id
            transcript_uri = "enst:" + transcript_stable_id

            self.triple(ordered_exon_uri, "a", "terms:EnsemblOrderedExon")
            self.triple(ordered_exon_uri, "a", "sio:SIO_001261")
            self.triple(ordered_exon_uri, "sio:SIO_000628", exon_uri)
            self.triple(ordered_exon_uri, "sio:SIO_000300", rank)

            self.triple(transcript_uri, "so:has_part", exon_uri)
            self.triple(transcript_uri, "sio:SIO_000974", ordered_exon_uri)
        self.output_file = sys.stdout
        f.close()
        return

    def rdfize_xref(self):
        gene = self.dbs["gene"]
        transcript = self.dbs["transcript"]
        translation = self.dbs["translation"]
        xref = self.dbs["xref"]
        object_xref = self.dbs["object_xref"]
        external_db = self.dbs["external_db"]
        f = open("xref.ttl", mode="w")
        self.output_file = f
        self.output_prefixes()
        for id in object_xref:
            xref_id = object_xref[id][2]
            subject_id = object_xref[id][0]
            subject_type = object_xref[id][1]
            if subject_type == "Gene":
                subject_url = "ensg:" + gene[subject_id][6]
            elif subject_type == "Transcript":
                subject_url = "enst:" + transcript[subject_id][7]
            elif subject_type == "Translation":
                subject_url = "ensp:" + translation[subject_id][1]
            else:
                continue
            # xref_node = Bnode()
            # xref_node.add(("terms:id_of", quote(external_db[xref[xref_id][0]][0])))  # FIXME
            # xref_node.add(("dcterms:identifier", quote(xref[xref_id][1])))
            # self.triple(subject_url, "rdfs:seeAlso", xref_node.serialize())
            external_db_id = xref[xref_id][0]
            external_db_code = external_db[external_db_id][0]
            if self.xref_url_dic.get(external_db_id, "") != "":
                xref_url = self.xref_url_dic[external_db_id] + xref[xref_id][1]
                self.triple(subject_url, "rdfs:seeAlso", "<"+xref_url+">")
                if external_db_code not in self.xrefed_dbs[subject_type]:
                    self.xrefed_dbs[subject_type][external_db_code] = [subject_url, xref_url, 0]
                self.xrefed_dbs[subject_type][external_db_code][2] += 1
            else:
                if external_db_code not in self.not_xrefed_dbs[subject_type]:
                    self.not_xrefed_dbs[subject_type][external_db_code] = [subject_url, xref[xref_id][1], 0]
                self.not_xrefed_dbs[subject_type][external_db_code][2] += 1
        self.output_file = sys.stdout
        f.close()
        return

    def output_prefixes(self):
        for prefix in Ensembl2turtle.prefixes:
            self.triple("@prefix", prefix[0], prefix[1])
        return

    def output_turtle(self):
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

        xref_report = {"xref": self.xrefed_dbs, "not_xref": self.not_xrefed_dbs}
        with open("xref_report.json", "w") as f:
            json.dump(xref_report, f, indent=4)

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
