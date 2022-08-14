# *_* coding: utf-8 *_*
# @Time     : 10/15/2021 11:16 AM
# @Author   : zong hui
# @object   : some useful functions related to PubTator.gov database


import os
import sys
import csv
import json
import getopt
import pandas as pd
from collections import defaultdict
from RecognizeGO import exact_string_match


def GOIDNameConverter(GO_info_file = '../knol/GO/go.info.tsv'):
    GO_ID2Name = dict()
    GO_Name2ID = dict()
    GO_IDName2ID = dict()
    df = pd.read_csv(GO_info_file, sep='\t')
    for GO_ID, GO_Name in zip(df['id'], df['name']):
        GO_ID2Name[GO_ID] = GO_Name
        GO_Name2ID[GO_Name] = GO_ID
        GO_IDName2ID[GO_ID] = GO_ID
        GO_IDName2ID[GO_Name] = GO_ID
    return GO_IDName2ID

def get_abstracts(pubtator_biocjson_path):
    pmid2abstract = dict()
    with open(pubtator_biocjson_path, 'r') as f:
        for line in f:
            data = json.loads(line)
            pmid = data['pmid']
            for passage in data['passages']:
                if passage['infons']['type'] == 'abstract':
                    abstract_text = passage.get('text', '')
            pmid2abstract[pmid] = abstract_text
    return pmid2abstract

def get_sentences(pubtator_sent_path):
    sentid2sentence = dict()
    with open(pubtator_sent_path, 'r') as f:
        f.readline()
        for line in f:
            l = line.strip('\n').split('\t')
            sentid2sentence[l[1]] = l[2]
    return sentid2sentence

def pubtator2go(case, pubtator_biocjson_path=None, pubtator_sent_path=None, pubmed2go_path=None):
    GO_IDName2ID = GOIDNameConverter()
    terms = list(GO_IDName2ID.keys())

    if pubtator_biocjson_path:
        # pubmed2abstract
        id2text = get_abstracts(pubtator_biocjson_path)
    elif pubtator_sent_path:
        # sentid2sentence
        id2text = get_sentences(pubtator_sent_path)
    ids = list(id2text.keys())
    texts = list(id2text.values())
    matched_terms = exact_string_match(terms, texts)
    ids2go = defaultdict(list)
    for id_, go_term in zip(ids, matched_terms):
        ids2go[id_] = [GO_IDName2ID[g] for g in go_term]

    print("[case]:", case, "[texts]:", len(texts))
    with open(pubmed2go_path, 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["PubMed ID", "GO Terms"])
        for id_, goid in ids2go.items():
            writer.writerow([id_, ";".join(goid)])
    return "done!"



def main(argv):
    # parse parameters
    case = ''
    try:
        opts, args = getopt.getopt(argv, "hc:", ["case="])
    except getopt.GetoptError:
        print('python pubtator2GO.py -c <case>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python pubtator2GO.py -c <case>')
            sys.exit()
        elif opt in ("-c", "--case"):
            case = arg
    
    # recognize gene ontology
    pubtator_sent_path = f'../case/{case}/PubTator/{case}.pubtator.bioclist.sent.txt'
    save_folder = f'../case/{case}/BiologicalEntity/'
    if not os.path.exists(save_folder):
        os.mkdir(save_folder)
    pubmed2go_path = f'../case/{case}/BiologicalEntity/{case}_pubmed2go.csv'
    pubtator2go(case, pubtator_sent_path=pubtator_sent_path, pubmed2go_path=pubmed2go_path)

if __name__ == "__main__":
    main(sys.argv[1:])

