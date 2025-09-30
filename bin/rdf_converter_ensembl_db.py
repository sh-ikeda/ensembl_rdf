import gzip
import sys
import os
import psutil
import json
import re
import datetime
import urllib.parse

input_dir = "./"


def quote(string):
    # Enclose a string with double-quotation marks.
    # Double-quotation marks within the input string are escaped with backslashes.
    return "\"" + string.replace("\"", "\\\"") + "\""


def escape(string):
    #return re.sub(r"([()])", r"\\\1", string)
    return urllib.parse.quote(string)


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
        ['ensgloss:', '<http://ensembl.org/glossary/>'],
        ['terms:', '<http://rdf.ebi.ac.uk/terms/ensembl/>'],
        ['ense:', '<http://rdf.ebi.ac.uk/resource/ensembl.exon/>'],
        ['ensp:', '<http://rdf.ebi.ac.uk/resource/ensembl.protein/>'],
        ['enst:', '<http://rdf.ebi.ac.uk/resource/ensembl.transcript/>'],
        ['ensi:', '<http://identifiers.org/ensembl/>'],
        ['taxonomy:', '<http://identifiers.org/taxonomy/>'],
        ['uniprot:', '<http://purl.uniprot.org/uniprot/>'],
        ['refseq:', '<http://identifiers.org/refseq/>'],
        ['sio:', '<http://semanticscience.org/resource/>'],
        ['skos:', '<http://www.w3.org/2004/02/skos/core#>']
    ]

    # Keys are `code` of the attrib_type table to be used as transcript flags
    transcript_flags = {
        "gencode_basic": {
            "GENCODE basic": "ensgloss:ENSGLOSSARY_0000020",
            "1": "ensgloss:ENSGLOSSARY_0000020"
        },
        "appris": {
            "principal1": "ensgloss:ENSGLOSSARY_0000013",
            "principal2": "ensgloss:ENSGLOSSARY_0000014",
            "principal3": "ensgloss:ENSGLOSSARY_0000015",
            "principal4": "ensgloss:ENSGLOSSARY_0000016",
            "principal5": "ensgloss:ENSGLOSSARY_0000017",
            "alternative1": "ensgloss:ENSGLOSSARY_0000018",
            "alternative2": "ensgloss:ENSGLOSSARY_0000019"
        },
        "TSL": {
            "tsl1": "ensgloss:ENSGLOSSARY_0000006",
            "tsl2": "ensgloss:ENSGLOSSARY_0000007",
            "tsl3": "ensgloss:ENSGLOSSARY_0000008",
            "tsl4": "ensgloss:ENSGLOSSARY_0000009",
            "tsl5": "ensgloss:ENSGLOSSARY_0000010",
            "tslNA": "ensgloss:ENSGLOSSARY_0000011"
        },
        "is_canonical": {"1": "ensgloss:ENSGLOSSARY_0000023"},
        # "remark": {"MANE_select": "ensgloss:ENSGLOSSARY_0000365"},
        "MANE_Select": {},
        "MANE_Plus_Clinical": {},
        "mRNA_start_NF": {"1": "ensgloss:ENSGLOSSARY_0000021"},
        "mRNA_end_NF": {"1": "ensgloss:ENSGLOSSARY_0000022"},
        "cds_start_NF": {"1": "ensgloss:ENSGLOSSARY_0000021"},
        "cds_end_NF": {"1": "ensgloss:ENSGLOSSARY_0000022"}
    }

    hco_chr_names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                     "11", "12", "13", "14", "15", "16", "17", "18",
                     "19", "20", "21", "22", "X", "Y", "MT"]

    base_dir = os.path.dirname(os.path.abspath(__file__)) + "/../"

    def __init__(self, input_dbinfo_file):
        self.dbinfo = self.load_dbinfo(input_dbinfo_file)
        self.dbs = self.load_dbs()
        # self.taxonomy_id = self.get_taxonomy_id()
        self.ensembl_version = self.get_ensembl_version()
        # self.production_name = self.get_production_name()
        self.species_id2taxonomy_id = self.get_species_id2taxonomy_id()
        # print(self.species_id2taxonomy_id)
        self.species_id2production_name = self.get_species_id2production_name()
        self.xref_url_dic = {}
        self.xref_prefix_dic = {}
        self.init_xref_url_dic()
        self.xrefed_dbs = {"Gene": {}, "Transcript": {}, "Translation": {}}
        self.not_xrefed_dbs = {"Gene": {}, "Transcript": {}, "Translation": {}}
        self.output_file = sys.stdout
        self.biotype_url_dic = {}
        self.init_biotype_url_dic()

    def init_biotype_url_dic(self):
        biotype_url_dic_tsv = "ontology/biotype_url.tsv"
        with open(Ensembl2turtle.base_dir+biotype_url_dic_tsv, "r") as input_table:
            line = input_table.readline()
            while (line):
                line = line.rstrip('\n')
                sep_line = line.split('\t')
                if sep_line[1] != "":
                    self.biotype_url_dic[sep_line[0]] = sep_line[1]
                line = input_table.readline()
        return

    def init_xref_url_dic(self):
        xref_url_dic_tsv = "config/external_db_url.tsv"

        with open(Ensembl2turtle.base_dir+xref_url_dic_tsv, "r") as input_table:
            line = input_table.readline()
            while (line):
                line = line.rstrip('\n')
                sep_line = line.split('\t')
                if sep_line[1] != "":
                    self.xref_url_dic[sep_line[0]] = sep_line[1]
                if sep_line[2] != "":
                    self.xref_prefix_dic[sep_line[0]] = sep_line[2]
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

    def get_species_id2production_name(self):
        return {v[0]: v[2] for v in self.dbs["meta"].values() if v[1] == 'species.production_name'}

    def get_species_id2taxonomy_id(self):
        return {v[0]: v[2] for v in self.dbs["meta"].values() if v[1] == 'species.taxonomy_id'}

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
            sbj = "ensg:" + escape(gene[id][6])
            xref_id = gene[id][4]
            seq_region_id = gene[id][7]

            self.triple(sbj, "a", "terms:EnsemblGene")
            biotype = gene[id][0]
            if biotype not in self.biotype_url_dic:
                print(f'Warning: Unknown biotype `{biotype}`', file=sys.stderr)
            else:
                self.triple(sbj, "a", self.biotype_url_dic[biotype])
                self.triple(sbj, "terms:has_biotype", self.biotype_url_dic[biotype])
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
            self.triple(sbj, "obo:RO_0002162", "taxonomy:"+self.seq_region_id_to_taxonomy_id(seq_region_id))

            # synonym
            synonyms = ", ".join([quote(v[0]) for v in external_synonym.get(xref_id, [])])
            if len(synonyms) >= 1:
                self.triple(sbj, "skos:altLabel", synonyms)

            # location
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
        unknown_flags = set()
        f = open("transcript.ttl", mode="w")
        self.output_file = f
        self.output_prefixes()
        for id in transcript:
            stable_id = transcript[id][7]
            sbj = "enst:" + escape(stable_id)
            xref_id = transcript[id][4]

            self.triple(sbj, "a", "terms:EnsemblTranscript")
            biotype = transcript[id][5]
            if biotype not in self.biotype_url_dic:
                print(f'Warning: Unknown biotype `{biotype}`', file=sys.stderr)
            else:
                self.triple(sbj, "a", self.biotype_url_dic[biotype])
                self.triple(sbj, "terms:has_biotype", self.biotype_url_dic[biotype])
            if xref_id == "\\N":
                label = stable_id  # Substitute ID for label
            else:
                label = xref[xref_id][2]
            self.triple(sbj, "rdfs:label", quote(label))
            self.triple(sbj, "dcterms:identifier", quote(stable_id))
            self.triple(sbj, "so:transcribed_from", "ensg:"+escape(gene[transcript[id][0]][6]))
            translates_to = transcript[id][6]
            if translates_to != "\\N":
                self.triple(sbj, "so:translates_to", "ensp:"+escape(translation[translates_to][1]))

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
                            statement = "<http://rdf.ebi.ac.uk/resource/ensembl.transcript/#_" + stable_id + "-has_transcript_flag-"+attrib_val+">"
                            self.triple(statement, "a", "rdf:Statement")
                            self.triple(statement, "rdf:subject", sbj)
                            self.triple(statement, "rdf:predicate", "terms:has_transcript_flag")
                            self.triple(statement, "rdf:object", flag_dic[attrib_code][attrib_val])
                            self.triple(statement, "rdfs:comment", quote(comment))
                    # elif attrib_code == "remark":
                    #     if attrib_val != "MANE_select":
                        #     continue
                        # else:
                    ## For the MANE transcripts, triples with a versioned transcript is added.
                    elif attrib_code == "MANE_Select" or attrib_code == "MANE_Plus_Clinical":
                        if attrib_code == "MANE_Select":
                            ensgloss_term = "ensgloss:ENSGLOSSARY_0000365"
                        else:
                            ensgloss_term = "ensgloss:ENSGLOSSARY_0000375"
                        version = transcript[id][9]
                        versioned_id = escape(stable_id) + "." + version
                        versioned_sbj = "enst:" + versioned_id
                        self.triple(sbj, "terms:has_transcript_flag", ensgloss_term)
                        self.triple(sbj, "terms:has_versioned_transcript", versioned_sbj)
                        self.triple(versioned_sbj, "a", "terms:VersionedTranscript")
                        self.triple(versioned_sbj, "terms:has_version", version)
                        self.triple(versioned_sbj, "dcterms:identifier", quote(versioned_id))
                        self.triple(versioned_sbj, "terms:has_transcript_flag", ensgloss_term)
                        counterpart = "refseq:" + attrib_val
                        self.triple(versioned_sbj, "terms:has_counterpart", counterpart)
                        self.triple(re.sub(r"\.[0-9]+$", "", counterpart), "terms:has_versioned_transcript", counterpart)
                        continue
                    # self.triple(sbj, "terms:has_transcript_flag", flag_dic[attrib_code][attrib_val])
                    try:
                        self.triple(sbj, "terms:has_transcript_flag", flag_dic[attrib_code][attrib_val])
                    except KeyError as e:
                        print(f"Warning: KeyError: {e}; {sbj} {attrib_code} {attrib_val}", file=sys.stderr)
                else:
                    if attrib[0] not in unknown_flags:
                        print(f"Warning: Attribute not treated as flag: {attrib[0]} {attrib_code}", file=sys.stderr)
                        unknown_flags.add(attrib[0])
        self.output_file = sys.stdout
        f.close()
        return

    def seq_region_id_to_taxonomy_id(self, seq_region_id):
        seq_region = self.dbs["seq_region"]
        coord_system = self.dbs["coord_system"]
        coord_system_id = seq_region[seq_region_id][1]
        species_id = coord_system[coord_system_id][0]
        return self.species_id2taxonomy_id[species_id]

    def seq_region_id_to_production_name(self, seq_region_id):
        seq_region = self.dbs["seq_region"]
        coord_system = self.dbs["coord_system"]
        coord_system_id = seq_region[seq_region_id][1]
        species_id = coord_system[coord_system_id][0]
        return self.species_id2production_name[species_id]

    def seq_region_id_to_chr(self, seq_region_id):
        seq_region = self.dbs["seq_region"]
        coord_system = self.dbs["coord_system"]
        chromosome_name = seq_region[seq_region_id][0]
        coord_system_id = seq_region[seq_region_id][1]
        production_name = self.seq_region_id_to_production_name(seq_region_id)
        chromosome_urls = []
        # e.g. "GRCm38"
        coord_system_version = coord_system[coord_system_id][2]
        # e.g. <http://rdf.ebi.ac.uk/resource/ensembl/109/mus_musculus/GRCm38/Y>
        chromosome_url = "<http://rdf.ebi.ac.uk/resource/ensembl/"+self.ensembl_version+"/"+production_name+"/"+coord_system_version+"/"+chromosome_name+">"
        # For LRG, <http://rdf.ebi.ac.uk/resource/ensembl/109/homo_sapiens/LRG_1>">"
        if coord_system[coord_system_id][1] == "lrg":
            chromosome_url = "<http://rdf.ebi.ac.uk/resource/ensembl/"+self.ensembl_version+"/"+production_name+"/"+chromosome_name+">"
        chromosome_urls.append(chromosome_url)

        if self.seq_region_id_to_taxonomy_id(seq_region_id) == "9606":
            if chromosome_name in Ensembl2turtle.hco_chr_names:
                hco_url = "<http://identifiers.org/hco/"+chromosome_name+"/"+coord_system_version+">"
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
            sbj = "ensp:" + escape(translation[id][1])

            self.triple(sbj, "a", "terms:EnsemblProtein")
            self.triple(sbj, "dcterms:identifier", quote(translation[id][1]))
            self.triple(sbj, "so:translation_of", "enst:"+escape(transcript[translation[id][0]][7]))
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
            sbj = "ense:" + escape(exon[id][3])

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
            ordered_exon_uri = "<http://rdf.ebi.ac.uk/resource/ensembl.transcript/"+escape(transcript_stable_id)+"#Exon_"+rank+">"
            exon_uri = "ense:" + escape(exon_stable_id)
            transcript_uri = "enst:" + escape(transcript_stable_id)

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
                subject_url = "ensg:" + escape(gene[subject_id][6])
            elif subject_type == "Transcript":
                subject_url = "enst:" + escape(transcript[subject_id][7])
            elif subject_type == "Translation":
                subject_url = "ensp:" + escape(translation[subject_id][1])
            else:
                continue
            # xref_node = Bnode()
            # xref_node.add(("terms:id_of", quote(external_db[xref[xref_id][0]][0])))  # FIXME
            # xref_node.add(("dcterms:identifier", quote(xref[xref_id][1])))
            # self.triple(subject_url, "rdfs:seeAlso", xref_node.serialize())
            external_db_id = xref[xref_id][0]
            external_db_code = external_db[external_db_id][0]
            if self.xref_url_dic.get(external_db_id, "") != "":
                dbprimary_acc = xref[xref_id][1]
                if external_db_id in self.xref_prefix_dic:
                    dbprimary_acc = dbprimary_acc.replace(self.xref_prefix_dic[external_db_id], "")
                xref_url = self.xref_url_dic[external_db_id] + urllib.parse.quote(dbprimary_acc, safe=":/")
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

        with open("xref_report.tsv", "w") as f:
            cwd = os.getcwd()
            dir_prod_name = re.sub(r"(.*/)|(_core_[^/]+$)", "", cwd)
            for subject_type, dbs in self.xrefed_dbs.items():
                for db in dbs:
                    print(dir_prod_name, "xref", subject_type, db,
                          dbs[db][0], dbs[db][1], dbs[db][2], sep="\t", file=f)
            for subject_type, dbs in self.not_xrefed_dbs.items():
                for db in dbs:
                    print(dir_prod_name, "unknown", subject_type, db,
                          dbs[db][0], dbs[db][1], dbs[db][2], sep="\t", file=f)

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
