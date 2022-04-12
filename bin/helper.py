# *_* coding: utf-8 *_*
# @Time     : 12/12/2020 12:20 PM
# @Author   : zong hui
# @object   : load some files


import csv
import gzip
import pandas as pd


def get_GeneIDNameMapping(infile="../knol/HumanGeneInformation/HumanGeneInformation.txt"):
    """get Gene ID-Name Mapping from downloaded files from HumanGeneInformation_file."""
    Name2IDs = dict()
    with open(infile, "r") as f:
        f.readline()
        for line in f:
            l = line.strip().split("\t")
            id_, name = l[0], l[1]
            if name not in Name2IDs:
                Name2IDs[name] = list()
            Name2IDs[name].append(id_)

    Gene_id2name = dict()
    Gene_name2id = dict()
    Gene_altid2id = dict()
    for name, ids in Name2IDs.items():
        if len(ids) == 1:
            Gene_name2id[name] = ids[0]
            Gene_id2name[ids[0]] = name
            Gene_altid2id[ids[0]] = ids[0]
        elif len(ids) > 1:
            ids_sorted = sorted(ids, key=lambda x: int(x), reverse=True)
            for i in range(len(ids_sorted)):
                Gene_name2id[name] = ids_sorted[0]
                Gene_id2name[ids_sorted[i]] = name
                Gene_altid2id[ids_sorted[i]] = ids_sorted[0]
    return Gene_id2name, Gene_name2id, Gene_altid2id


def get_GOIDNameMapping(infile="../knol/GO/go.info.csv"):
    GOinfo = pd.read_csv(infile, keep_default_na=False)
    GO_id2name = dict(zip(GOinfo.id, GOinfo.name))
    GO_name2id = dict(zip(GOinfo.name, GOinfo.id))
    GO_id2level = dict(zip(GOinfo.id, GOinfo.level))
    GO_id2namespace = dict(zip(GOinfo.id, GOinfo.namespace))
    GO_id2children = {}
    for id_, children_ in zip(GOinfo.id, GOinfo.children):
        GO_id2children[id_] = children_.split(';')
    GO_altid2id = {}
    for go_id, alt_ids in zip(GOinfo.id, GOinfo.alt_ids):
        if alt_ids:
            alt_ids_list = alt_ids.split(';')
            for alt_id in alt_ids_list:
                GO_altid2id[alt_id] = go_id
    return GO_id2name, GO_name2id, GO_id2level, GO_id2children, GO_id2namespace, GO_altid2id


def get_NCBI_GeneGO2PMID(infile="../knol/NCBI_EntrezGene/gene2go.gz"):
    # GO annotations from the gene2go.gz nfile avalable at NCBI EntrezGene ftp
    NCBI_GeneGO2PMID = dict()
    with gzip.open(infile, 'rb') as f:
        f.readline()
        for line in f:
            l = line.decode().split("\t")
            gene_id, go_id, pubmed = l[1], l[2], l[6]
            if pubmed != '-':
                NCBI_GeneGO2PMID[(gene_id, go_id)] = pubmed.split('|')
    return NCBI_GeneGO2PMID


def get_MSigDB_GO2Gene(MSigDB_file="../knol/MSigDB/c5.go.v7.2.symbols.gmt"):
    MSigDB_GO2Gene = {}
    with open(MSigDB_file, 'r') as f:
        for line in f:
            l = line.strip().split("\t")
            MSigDB_GO2Gene[l[0].replace("GO_","").replace("_"," ").lower()] = l[2:]
    return MSigDB_GO2Gene
