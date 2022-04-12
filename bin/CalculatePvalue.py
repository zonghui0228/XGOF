# *_* coding: utf-8 *_*
# @Time    : 2020/11/27 21:16
# @Author  : zong hui


# load packages
import scipy
import numpy as np
import pandas as pd
from tqdm import tqdm
import scipy.stats as stats
from collections import defaultdict

# load our modules
import config

basicConfig = config.basicConfig
from helper import get_GeneIDNameMapping, get_GOIDNameMapping
from helper import get_MSigDB_GO2Gene, get_NCBI_GeneGO2PMID

# Gene ID -- Gene Name
Gene_id2name, Gene_name2id, Gene_altid2id = get_GeneIDNameMapping(basicConfig.HumanGeneInformation)
# GO ID -- GO Name
GO_id2name, GO_name2id, GO_id2level, GO_id2children, GO_id2namespace, GO_altid2id = get_GOIDNameMapping(basicConfig.GO_INFO)
# MSigDB file
MSigDB_GO2Gene = get_MSigDB_GO2Gene(basicConfig.MSigDB_C5_GO_PATH)
# NCBI gene2go
NCBI_GeneGO2PMID = get_NCBI_GeneGO2PMID(basicConfig.NCBI_Entrez_Gene2GO)


class CalculatePvalue:
    """
    Calculate Pvalue, Adjusted Pvalue of Gene and GO based on fisher exact test.
    Args:
        GeneMapping (:obj: 'dict'):
            the Gene mapping result, the key is gene, and the value is related pubmed id.
        GoMapping (:obj: 'dict'):
            the GO mapping result, the key is GO, and the value is related pubmed id.
    """

    def __init__(self, mapping_gene2pmid, mapping_go2pmid):
        self.GeneMapping = self.get_GeneMapping(mapping_gene2pmid)  # input Gene mapping
        self.GoMapping = self.get_GoMapping(mapping_go2pmid)  # input GO mapping

        self.pmid = list(set([p for k, v in self.GeneMapping.items() for p in v]))  # get the initial PMID list.
        self.genes = list(self.GeneMapping.keys())  # get initial Genes
        self.gos = list(self.GoMapping.keys())  # get initial GOs

        self.GOF = defaultdict(dict)
        self.GOF_Genes = list()
        self.GOF_GOs = list()
        self.GeneGoEnrichment = pd.DataFrame(
            columns=["Letter", "Gene_ID", "Gene_Name", "GO_ID", "GO_Namespace", "GO_Name", "GOF",
                     "a", "b", "c", "d", "Pvalue", "Adjusted_Pvalue", "Enrichment_Score", "GO_Level", "GO_Children",
                     "GSEA_MSigDB", "NCBI_Entrez"])

        self.calculate_pvalue_adjustedpvalue()
        self.calculate_enrichment_score()
        self.add_evidence()
        self.filter_and_sort()
        

    def get_GeneMapping(self, mapping_gene2pmid):
        GeneMapping = {}
        with open(mapping_gene2pmid, 'r') as f:
            f.readline()
            for line in f:
                l = line.strip().split('\t')
                GeneMapping[l[0]] = l[2].split(';')
        return GeneMapping

    def get_GoMapping(self, mapping_go2pmid):
        GoMapping = {}
        with open(mapping_go2pmid, 'r') as f:
            f.readline()
            for line in f:
                l = line.strip().split('\t')
                GoMapping[l[0]] = l[2].split(';')
        return GoMapping

    def calculate_pvalue_adjustedpvalue(self):
        tqdm_genes = tqdm(self.genes, ncols=80)
        entity = self.GoMapping
        D = dict()
        for gene in tqdm_genes:
            tqdm_genes.set_description("Calculate Pvalue and Adjusted Pvalue")
            gene_pmid = set(self.GeneMapping.get(gene))

            gene_entity2pvalue = dict()
            gene_pvalues = list()
            for e in entity:
                e_pmid = set(entity.get(e))
                a = len(gene_pmid & e_pmid)
                b = len(gene_pmid - e_pmid)
                c = len(e_pmid - gene_pmid)
                d = len(self.pmid) - a - b - c
                if a>0:
                    oddsration, pvalue = stats.fisher_exact([[a, b], [c, d]], alternative='two-sided')
                    if pvalue < 1.0:
                        gene_pvalues.append(pvalue)
                        gene_entity2pvalue[e] = pvalue
                        D[(gene, e)] = {"Gene_ID": gene,
                                        "a": a,
                                        "b": b,
                                        "c": c,
                                        "d": d,
                                        "Pvalue": pvalue}
                        D[(gene, e)]["GO_ID"] = e

            for e, e_pvalue in gene_entity2pvalue.items():
                numerator, denominator = 0, 0
                for k, k_pvalue in gene_entity2pvalue.items():
                    if k_pvalue < 1.0:
                        numerator += 1
                    if k_pvalue < 1.0 and k_pvalue >= e_pvalue:
                        denominator += 1
                adjusted_pvalue = min(1.0, e_pvalue * float(numerator) / float(denominator))
                D[(gene, e)]["Adjusted_Pvalue"] = adjusted_pvalue

        tqdm_genes.close()
        self.GeneGoEnrichment = self.GeneGoEnrichment.append(list(D.values()), ignore_index=True)
        D = dict()

    def calculate_enrichment_score(self):
        adjusted_pvalue = self.GeneGoEnrichment['Adjusted_Pvalue'].tolist()
        adjusted_pvalue_minmum = min([i for i in adjusted_pvalue if i>0])
        enrichment_score = [-np.log(p) if p!=0 else -np.log(adjusted_pvalue_minmum)for p in adjusted_pvalue]
        self.GeneGoEnrichment["Enrichment_Score"] = enrichment_score
        self.Pvalue2AP =  pd.DataFrame(self.GeneGoEnrichment, columns=['Gene_ID', 'GO_ID', 'Pvalue', 'Adjusted_Pvalue'])

    def save_Pvalue2AdjustedPvalue(self, p2ap_save_path=None):
        if p2ap_save_path:
            self.Pvalue2AP.to_csv(p2ap_save_path, index=False)
        return 'done!'

    def add_evidence(self):
        # GOF: Letter, Gene Name, GO Name, GO level, GO namespace, GO children
        self.GeneGoEnrichment["Gene_Name"] = self.GeneGoEnrichment["Gene_ID"].apply(lambda x: Gene_id2name.get(x, x))
        self.GeneGoEnrichment["Letter"] = self.GeneGoEnrichment["Gene_Name"].apply(lambda x: x.upper()[0])
        self.GeneGoEnrichment["GO_Name"] = self.GeneGoEnrichment["GO_ID"].apply(lambda x: GO_id2name.get(x, x))
        self.GeneGoEnrichment["GO_Level"] = self.GeneGoEnrichment["GO_ID"].apply(lambda x: GO_id2level.get(x, x))
        self.GeneGoEnrichment["GO_Children"] = self.GeneGoEnrichment["GO_ID"].apply(lambda x: len(GO_id2children.get(x, [])))
        self.GeneGoEnrichment["GO_Namespace"] = self.GeneGoEnrichment["GO_ID"].apply(lambda x: GO_id2namespace.get(x, x))
        # GeneGOname = ["{}[{}]".format(gene_name, go_name) for (gene_name, go_name) in zip(self.GeneGoEnrichment["Gene_Name"], self.GeneGoEnrichment["GO_Name"])]
        # self.GeneGoEnrichment["GOF"] = np.array(GeneGOname)

        # GSEA
        MSigDB_evidence = [1 if gene_name in MSigDB_GO2Gene.get(go_name.lower(),[]) else 0 for (go_name, gene_name) in
                         zip(self.GeneGoEnrichment["GO_Name"], self.GeneGoEnrichment["Gene_Name"])]
        self.GeneGoEnrichment["GSEA_MSigDB"] = np.array(MSigDB_evidence)

        # PubMed
        NCBI_evidence = [len(NCBI_GeneGO2PMID.get((gene_id, go_id),[])) for (gene_id, go_id) in zip(self.GeneGoEnrichment["Gene_ID"], self.GeneGoEnrichment["GO_ID"])]
        self.GeneGoEnrichment["NCBI_Entrez"] = np.array(NCBI_evidence)

    def filter_and_sort(self):
        # 1. filter the Gene-GO, df[(df['a']>=5) | (df['GSEA_MSigDB']!=0) | (df['NCBI_Entrez']!=0)]
        self.GeneGoEnrichment = self.GeneGoEnrichment.loc[self.GeneGoEnrichment["Adjusted_Pvalue"] <= 0.05]
        self.GeneGoEnrichment = self.GeneGoEnrichment[(self.GeneGoEnrichment['a']>=5) | (self.GeneGoEnrichment['GSEA_MSigDB']!=0) | (self.GeneGoEnrichment['NCBI_Entrez']!=0)]

        # 2. sorted in descending order by Adjusted_Pvalue.
        self.GeneGoEnrichment = self.GeneGoEnrichment.sort_values('Adjusted_Pvalue', ascending=True)

        # 3. obtain GOF
        for (gene, go, ap) in zip(self.GeneGoEnrichment["Gene_ID"], self.GeneGoEnrichment["GO_ID"],self.GeneGoEnrichment["Adjusted_Pvalue"]):
            self.GOF[gene][go] = ap
        self.GOF_Genes = list(set(self.GeneGoEnrichment["Gene_ID"].tolist()))
        self.GOF_GOs = list(set(self.GeneGoEnrichment["GO_ID"].tolist()))

    def get_most_related_term(self, term="BCL2", term_type="Gene", max_number=10):
        res = []
        if term_type == "Gene":
            res = self.GeneGoEnrichment.loc[self.GeneGoEnrichment["Gene_Name"] == term].sort_values('Adjusted_Pvalue',ascending=False)[:max_number]
        if term_type == "GO":
            res = self.GeneGoEnrichment.loc[self.GeneGoEnrichment["GO_ID"] == term].sort_values('Adjusted_Pvalue',ascending=False)[:max_number]
        return res

    def save_GOF(self, GOF_save_path=None):
        if GOF_save_path:
            self.GeneGoEnrichment.to_csv(GOF_save_path, index=False)
