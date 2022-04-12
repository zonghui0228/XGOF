# *_* coding: utf-8 *_*
# @Time     : 10/15/2021 11:16 AM
# @Author   : zong hui
# @object   : some useful functions related to PubTator.gov database


import csv
import json
from tqdm import tqdm
from collections import defaultdict



def pubtator2gene(case, pubtator_path, pubmed2gene_path):
    num_lines = sum(1 for line in open(pubtator_path, 'r'))
    print("case: {}, pmid: {}".format(case, num_lines))
    # 获得PMID2Gene
    pmid2gene = defaultdict(list)
    with open(pubtator_path, 'r') as f:
        for line in tqdm(f, total=num_lines, ncols=80):
            data = json.loads(line)
            pmid = data['pmid']
            genes = list()
            for passage in data['passages']:
                if passage['infons']['type'] == 'abstract':
                    # 摘要文本
                    abstract_text = passage.get('text', '')
                    # 注释信息
                    annotations = passage.get('annotations', [])
                    if annotations:
                        for annotation in annotations:
                            # 基因注释
                            if annotation['infons']['type'] == 'Gene':
                                genes.append(annotation['infons']['identifier'])
            pmid2gene[pmid] = list(set(genes))
    # 保存PMID2Gene
    with open(pubmed2gene_path, 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["PubMed ID", "Gene Terms"])
        for pmid, genes in pmid2gene.items():
            # 筛选保留人类基因
            human_gene = [gene for gene in genes if Gene_ID2Name.get(gene)]
            writer.writerow([pmid, ";".join(human_gene)])
    return "done!"


def pubtator2gene1(case, pubtator_anno_path, pubmed2gene_path):
    # 获取gene映射文件
    pmid, gene, sentenceid = set(), set(), set()
    sentid2gene = defaultdict(list)
    with open(pubtator_anno_path, 'r') as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')
            if l[4]=='Gene':
                pmid.add(l[0])
                gene.add(l[5])
                sentenceid.add(l[-1])
                sentid2gene[l[-1]].append(l[5])
    print("[case]: {}, [pmid]: {}, [sentence]: {}, [gene]: {}".format(case, len(pmid), len(sentenceid), len(gene)))
    # 保存PMID2Gene
    with open(pubmed2gene_path, 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["PubMed ID", "Gene Terms"])
        for sid, geneids in sentid2gene.items():
            # 筛选保留人类基因
            writer.writerow([sid, ";".join(list(set(geneids)))])
    return "done!"


if __name__ == "__main__":
    # ==================================================
    # cases = ['AAV2', 'COVID19', 'EBV', 'HBV', 'HIV', 'HPV', 'HTLV1', 'MCV', 'XMRV']
    # for case in cases:
    #     pubtator_anno_path = '../data/PubTator/virus/{}/{}.pubtator.bioclist.anno.txt'.format(case, case)
    #     pubmed2gene_path = '../data/BiologicalEntity/Gene2/virus/{}_pubmed2gene.csv'.format(case)
    #     pubtator2gene1(case, pubtator_anno_path, pubmed2gene_path)
        
    # ==================================================
    # cases = ['CHOL', 'COADREAD', 'ESCA', 'LIHC', 'PAAD', 'STAD']
    # for case in cases:
    #     pubtator_anno_path = '../data/PubTator/GIDB/{}/{}.pubtator.bioclist.anno.txt'.format(case, case)
    #     pubmed2gene_path = '../data/BiologicalEntity/Gene2/GIDB/{}_pubmed2gene.csv'.format(case)
    #     pubtator2gene1(case, pubtator_anno_path, pubmed2gene_path)
    
    # ==================================================
    # cases = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH',
    #          'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG',
    #          'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
    # for case in cases:
    #     pubtator_anno_path = '../data/PubTator/TCGA/{}/{}.pubtator.bioclist.anno.txt'.format(case, case)
    #     pubmed2gene_path = '../data/BiologicalEntity/Gene2/TCGA/{}_pubmed2gene.csv'.format(case)
    #     pubtator2gene1(case, pubtator_anno_path, pubmed2gene_path)

    print('done!')
