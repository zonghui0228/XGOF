# *_* coding: utf-8 *_*
# @Time     : 11/11/2020 12:59 PM
# @Author   : zong hui
# @object   : given a list of pmid, get the mapping of pmid to gene


import os
import csv
import gzip
from collections import defaultdict

from helper import getGeneIDNameMapping
Gene_ID2Name, Gene_Name2ID, Gene_altID2ID = getGeneIDNameMapping("../knol/HumanGeneInformation/HumanGeneInformation.txt")

def getGenefromPubtator(pmid_path, pubtator_path, pmid2gene_path):
    # 获取pmid列表
    with open(pmid_path, 'r') as f:
        pmid = [line.strip() for line in f]
    
    # 如果已经提取了gene，我们从相应的pmid.txt中读取
    if os.path.exists(pmid2gene_path):
        genes = list()
        with open(pmid2gene_path, 'r') as f:
            rows = list(csv.reader(f))
            for row in rows[1:]:
                genes.extend(row[1].split(';'))
        genes = list(set(genes))
        genes.remove('')
        
    # 如果还未提取gene，那开始提取
    else:
        pubmed2gene = defaultdict(list)
        with gzip.open(pubtator_path, 'r') as f:
            for line in f:
                l = line.decode().split('\t')
                if Gene_ID2Name.get(l[2]):
                    pubmed2gene[l[0]].append(l[2])

        genes = list()
        with open(pmid2gene_path, 'w', newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["PubMed ID", "Gene Terms"])
            for p in pmid:
                writer.writerow([p, ";".join(pubmed2gene.get(p,''))])
                genes.extend(pubmed2gene.get(p, []))
        genes = list(set(genes))
    print(len(pmid), len(genes))
    
    
if __name__ == "__main__":
    pmid_path = "../case/case-virus/HIV/HIV.pmid.txt"
    pubtator_path = "../../rawdata/PubTator/gene2pubtatorcentral.gz"
    pmid2gene_path = "../db/Gene/HIV_pubmed2gene.csv"
    getGenefromPubtator(pmid_path, pubtator_path, pmid2gene_path)

