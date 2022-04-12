# *_* coding: utf-8 *_*
# @Time    : 2020/11/27 13:50
# @Author  : zong hui


# load packages
import csv
from collections import defaultdict


# load our modules
import config
basicConfig = config.basicConfig
from helper import get_GeneIDNameMapping, get_GOIDNameMapping
# Gene ID -- Gene Name
Gene_id2name, Gene_name2id, Gene_altid2id = get_GeneIDNameMapping(basicConfig.HumanGeneInformation)
# GO ID -- GO Name
GO_id2name, GO_name2id, GO_id2level, GO_id2children, GO_id2namespace, GO_altid2id = get_GOIDNameMapping(basicConfig.GO_INFO)


class EntityMapping:
    """
    This is the EntityMapping class to get the Gene-Pubmed mapping and GO-Pubmed mapping.
    For a given pmid list, we extract related entities (Gene, GO)e, and convert into the predefined format.

    Args:
        pmid_list_path (:obj: 'string'):
            a file of pubmed id.
        pubmed2gene_path (:obj: 'string'):
            a file path which records the Gene entity related to Pubmed id.
        pubmed2go_path (:obj: 'string'):
            a file path which records the GO entity related to Pubmed id.
    """
    def __init__(self, pmid_list_path, pubmed2gene_path=None, pubmed2go_path=None):
        self.pmid_list_path = pmid_list_path
        self.pubmed2gene_path = pubmed2gene_path
        self.pubmed2go_path = pubmed2go_path

        with open(self.pmid_list_path, "r") as f1:
            self.input_pmid_list = [line.strip() for line in f1]

        self.pmid_list = list()
        self.PMIDEntitesMapping = dict()
        self.get_PMIDEntitesMappingWithFiles()

        self.GeneMapping = defaultdict(list)
        self.get_GeneMapping()
        self.GoMapping = defaultdict(list)
        self.get_GoMapping()

    # 获得PMID对于的Gene和GO，gene必须全部是出现在我们HumanGeneInformation中的。
    def get_PMIDEntitesMappingWithFiles(self):
        PMIDEntitesMapping = dict()
        input_pmid_dict = {pmid: "input_pmid" for pmid in self.input_pmid_list}

        for pmid in self.input_pmid_list:
            PMIDEntitesMapping[pmid] = {"Gene": list(), "GO": list()}

        with open(self.pubmed2gene_path, "r") as f1:
            rows = list(csv.reader(f1))
            for row in rows:
                pmid, this_genes = row[0], row[1].split(";")
                if input_pmid_dict.get(pmid):
                    PMIDEntitesMapping[pmid]["Gene"] = [Gene_altid2id.get(g, g) for g in this_genes if Gene_altid2id.get(g)]

        with open(self.pubmed2go_path, "r") as f2:
            rows = list(csv.reader(f2))
            for row in rows:
                pmid, this_gos = row[0], row[1].split(";")
                if input_pmid_dict.get(pmid):
                    PMIDEntitesMapping[pmid]["GO"] = [GO_altid2id.get(g, g) for g in this_gos if GO_id2name.get(g)]

        self.PMIDEntitesMapping = {k:v for k,v in PMIDEntitesMapping.items() if v["Gene"] and v["GO"]}
        self.pmid_list = self.PMIDEntitesMapping.keys() # 这是基因都是HumanGeneInformation里的PMID

    def get_GeneMapping(self, method="Pubtator"):
        for pmid in self.pmid_list:
            terms = self.PMIDEntitesMapping[pmid]["Gene"]
            for term in terms:
                self.GeneMapping[term].append(pmid)

    def get_GoMapping(self, method="MetaMap"):
        for pmid in self.pmid_list:
            terms = self.PMIDEntitesMapping[pmid]["GO"]
            for term in terms:
                self.GoMapping[term].append(pmid)

    def get_terms(self, pmid="10023698", type="Gene"):
        terms = self.PMIDEntitesMapping[pmid][type]
        return terms

    def save_PMIDEntitesMapping(self, save_path=None):
        with open(save_path, "w") as f:
            f.write("PMID\tGene_ID\tGO_ID\n")
            for pmid in self.PMIDEntitesMapping:
                genes = ";".join(self.PMIDEntitesMapping[pmid]["Gene"])
                gos = ";".join(self.PMIDEntitesMapping[pmid]["GO"])
                f.write("{}\t{}\t{}\n".format(pmid, genes, gos))

    def save_GeneMapping(self, save_path=None):
        with open(save_path, "w") as f:
            f.write(("Gene_ID\tGene_Terms\tPMID\n"))
            for k,v in self.GeneMapping.items():
                f.write("{}\t{}\t{}\n".format(k, Gene_id2name.get(k,k), ";".join(v)))

    def save_GoMapping(self, save_path=None):
        with open(save_path, "w") as f:
            f.write(("GO_ID\tGO_Terms\tPMID\n"))
            for k,v in self.GoMapping.items():
                f.write("{}\t{}\t{}\n".format(k, GO_id2name.get(k,k), ";".join(v)))
