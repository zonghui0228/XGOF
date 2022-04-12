# *_* coding: utf-8 *_*
# @Time    : 2020/11/27 21:16
# @Author  : zong hui


# load packages
import math
import pandas as pd
from tqdm import tqdm
from collections import defaultdict


# load our modules
import config
basicConfig = config.basicConfig
from helper import get_GeneIDNameMapping
# Gene ID -- Gene Name
Gene_id2name, Gene_name2id, Gene_altid2id = get_GeneIDNameMapping(basicConfig.HumanGeneInformation)


class GeneSimilarity:
    def __init__(self, GOF_file):
        self.GOF = dict
        self.GOF_Genes = list()
        self.GOF_GOs = list()
        self.GGSS_Node = set()
        self.GGSS_Edge = set()

        self.GeneGeneSimilarity = pd.DataFrame(
            columns=["Gene1_ID", "Gene1_Name", "Gene2_ID", "Gene2_Name", "Similarity_Score","Similarity_Score(MinMaxScaler)", "Similarity_Score(StandardScaler)"])
        self.get_GOF(GOF_file)
        self.calculate_GeneGeneSimilarity()
        self.normalize_GeneGeneSimilarity()

    def get_GOF(self, GOF_file):
        self.GOF = defaultdict(dict)
        df = pd.read_csv(GOF_file, keep_default_na=False)
        for (gene, go, ap) in zip(df["Gene_ID"], df["GO_ID"], df["Adjusted_Pvalue"]):
            self.GOF[str(gene)][go] = ap
        self.GOF_Genes = [str(g) for g in list(set(df["Gene_ID"].tolist()))]
        self.GOF_GOs = list(set(df["GO_ID"].tolist()))

    def calculate_GeneGeneSimilarity(self):
        # calculate gene similarity based on modified inner product.
        D = list()
        tqdm_GOF_Genes = tqdm(self.GOF_Genes, ncols=80)
        for gene1 in tqdm_GOF_Genes:
            tqdm_GOF_Genes.set_description("Calculate Gene-Gene Similarity Score")
            gene1_index = self.GOF_Genes.index(gene1)
            gene1_go = self.GOF[gene1]
            for gene2 in self.GOF_Genes[gene1_index+1:]:
                gene2_go = self.GOF[gene2]
                both_gene_go = set(gene1_go) & set(gene2_go)
                only_gene1_go = set(gene1_go) - set(gene2_go)
                only_gene2_go = set(gene2_go) - set(gene1_go)
                score = 0.0
                if len(both_gene_go) != 0:
                    numerator, demominator = 0.0, 0.0
                    for g in both_gene_go:
                        if self.GOF[gene1][g] > 0 and self.GOF[gene2][g] > 0:
                            numerator += math.log(self.GOF[gene1][g]) * math.log(self.GOF[gene2][g])
                    demominator = max(1.0, 0.5 * (len(only_gene1_go) + len(only_gene2_go)))
                    score = numerator / demominator
                    D.append({"Gene1_ID":gene1, "Gene1_Name":Gene_id2name.get(gene1,gene1), "Gene2_ID":gene2, "Gene2_Name":Gene_id2name.get(gene2,gene2), "Similarity_Score":score})
                    self.GGSS_Node.update({gene1, gene2})
                    self.GGSS_Edge.add((gene1, gene2))
        self.GeneGeneSimilarity = self.GeneGeneSimilarity.append(D, ignore_index=True)
        tqdm_GOF_Genes.close()

    def normalize_GeneGeneSimilarity(self):
        # normalize scores by scaling each score to a given range.
        # MinMaxScaler
        max_ = self.GeneGeneSimilarity["Similarity_Score"].agg(max)
        min_ = self.GeneGeneSimilarity["Similarity_Score"].agg(min)
        self.GeneGeneSimilarity["Similarity_Score(MinMaxScaler)"] = self.GeneGeneSimilarity["Similarity_Score"].apply(lambda x: (x - min_) / (max_ - min_))

        # StandardScaler
        mean_ = self.GeneGeneSimilarity['Similarity_Score'].agg("mean")
        std_ = self.GeneGeneSimilarity["Similarity_Score"].agg("std")
        self.GeneGeneSimilarity["Similarity_Score(StandardScaler)"] = self.GeneGeneSimilarity["Similarity_Score"].apply(lambda x: (x - mean_) / std_)

        # Sorted in descending order by Similarity Score.
        self.GeneGeneSimilarity.sort_values(by="Similarity_Score", inplace=True, ascending=False)

    def get_most_related_gene(self, term="CDK2", type="Name", max_number=10):
        res = pd.DataFrame()
        if type=="Name":
            res = self.GeneGeneSimilarity.loc[(self.GeneGeneSimilarity["Gene1_Name"]==term) | (self.GeneGeneSimilarity["Gene1_Name"]==term)][:max_number]
        if type=="ID":
            res = self.GeneGeneSimilarity.loc[(self.GeneGeneSimilarity["Gene1_ID"] == term) | (self.GeneGeneSimilarity["Gene1_ID"] == term)][:max_number]
        return res

    def save_GGSS(self, save_path=None, percentage=1):
        # filter GGSS data.
        GGSS_length = self.GeneGeneSimilarity.shape[0]
        GGSS_data = self.GeneGeneSimilarity[:int(percentage*GGSS_length)]
        num_node = len(set(GGSS_data["Gene1_ID"].tolist()+GGSS_data["Gene2_ID"].tolist()))
        num_edge = GGSS_data.shape[0]

        # save file
        GGSS_data.to_csv(save_path, index=False)
        return num_node, num_edge
