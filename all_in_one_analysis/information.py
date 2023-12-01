from Bio import ExPASy, SwissProt
from bioservices import UniProt
import re

def get_gene(uniprot_id):
    handle = ExPASy.get_sprot_raw(uniprot_id)
    record = SwissProt.read(handle)
    gene_name = record.gene_name
    names = []
    for i in record.gene_name:
        names.append(i['Name'])
    return names

def get_goterms(uniprot_id):
    u = UniProt(verbose=False)
    id = [i.strip() for i in u.search("P01011", columns = "go_id", frmt="tsv" ).split("\n")[1].split(";")]
    bp = [i.strip() for i in u.search("P01011", columns = "go_p", frmt="tsv" ).split("\n")[1].split(";")]
    cc = [i.strip() for i in u.search("P01011", columns = "go_c", frmt="tsv" ).split("\n")[1].split(";")]
    mf = [i.strip() for i in u.search("P01011", columns = "go_f", frmt="tsv" ).split("\n")[1].split(";")]
    goterms = {"id": id, "biological process": bp, "cellular component": cc, "molecular function": mf}
    return goterms

def get_description(uniprot_id):
    handle = ExPASy.get_sprot_raw(uniprot_id)
    record = SwissProt.read(handle)
    description = record.description
    description = re.findall(r'Full=(.+?);', description)
    return description

def get_comment(uniprot_id):
    handle = ExPASy.get_sprot_raw(uniprot_id)
    record = SwissProt.read(handle)
    comments = {}
    for c in record.comments:
        _c = c.split(":")
        comments[_c[0]] = ':'.join(_c[1:]).strip()
    return comments