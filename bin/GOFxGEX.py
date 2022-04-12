# *_* coding: utf-8 *_*
# @Time     : 3/11/2021 12:37 PM
# @Author   : zong hui
# @object   : matrix GOFxSeq(GOxSample) = GOF(GOxGene) x Seq(GenexSample)


import csv
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import defaultdict


class GOFxGEX:
    def __init__(self, GOF_file, GEX_file):
        self.GOF_file = GOF_file
        self.GEX_file = GEX_file
        
        self.GOF_matrix = pd.DataFrame()
        self.GOF_GOs = list()
        self.GEX_matrix = pd.DataFrame()
        self.GEX_Samples = list()
        
        self.sharedGene = list()
        self.GOF_sharedGene = pd.DataFrame()
        self.GEX_sharedGene = pd.DataFrame()
        self.GOFxGEX_matrix = pd.DataFrame()

        tqdm_ = tqdm(['tqdm'], ncols=80)
        for i in tqdm_:
            tqdm_.set_description("Calculate GOFxGEX Matirx")
            self.get_GOF()
            self.get_GEX()
            self.get_GOFxGEX()

    def get_overlapGene(self):
        GOF_info = list(csv.reader(open(self.GOF_file, 'r')))
        GOF_Genes = set([info[2] for info in GOF_info[1:]])
        with open(self.GEX_file, 'r') as f:
            GEX_Genes = set([line.split('\t')[0].split('|')[0] for line in f])
        self.OverlapGenes = sorted(list(GOF_Genes & GEX_Genes))

    def get_GOF(self):
        """生成GOF矩阵，行为GO Name,列为overlapGene Name"""
        GOF = pd.read_csv(self.GOF_file)
        self.GOF_GOs = GOF['GO_Name'].unique().tolist()
        self.GOF_matrix = pd.DataFrame(index=GOF['GO_Name'].unique(), columns=GOF['Gene_Name'].unique())
        for go, gene, es in zip(GOF['GO_Name'], GOF['Gene_Name'], GOF['Enrichment_Score']):
            self.GOF_matrix.loc[go, gene] = es
        self.GOF_matrix = self.GOF_matrix.fillna(0)
        self.GOF_matrix.index.name = 'GO({})xGene({})'.format(len(self.GOF_matrix.index), len(self.GOF_matrix.columns))

    def get_GEX(self):
        """将表达数据转换为overlapGenes x Samples矩阵"""
        self.GEX_matrix = pd.read_csv(self.GEX_file, sep='\t', index_col=0)
        self.GEX_matrix.index = [i.split('|')[0] for i in self.GEX_matrix.index]
        self.GEX_matrix.columns = [i[:15].replace('-','.') for i in self.GEX_matrix.columns]
        self.GEX_matrix.index.name = 'Gene({})xSample({})'.format(len(self.GEX_matrix.index), len(self.GEX_matrix.columns))
        
        self.GEX_Samples = self.GEX_matrix.columns.unique().tolist()
        self.Samples_tumor = [s for s in self.GEX_Samples if s[13:15] in ["0"+str(i) for i in range(1,10)]]
        self.Samples_normal = [s for s in self.GEX_Samples if s[13:15] in [str(i) for i in range(10,20)]]
        

    def get_GOFxGEX(self):
        self.sharedGene = list(set(self.GOF_matrix.columns.tolist()) & set(self.GEX_matrix.index.tolist()))
        self.GOF_sharedGene = self.GOF_matrix.loc[:, self.sharedGene]
        self.GEX_sharedGene = self.GEX_matrix.loc[self.sharedGene, :]

        GOFxGEX_values = np.dot(self.GOF_sharedGene.values, self.GEX_sharedGene.values)
        self.GOFxGEX_matrix = pd.DataFrame(GOFxGEX_values, index=self.GOF_sharedGene.index, columns=self.GEX_sharedGene.columns)
        self.GOFxGEX_matrix.index.name = 'GO({})xSample({})'.format(len(self.GOFxGEX_matrix.index), len(self.GOFxGEX_matrix.columns))
        
    def save_GOF2matrix(self, save_path):
        self.GOF_matrix.to_csv(save_path, sep='\t')

    def save_GEX2matrix(self, save_path):
        self.GEX_matrix.to_csv(save_path, sep='\t')

    def save_GOFxGEX(self, save_path):
        self.GOFxGEX_matrix.to_csv(save_path, sep='\t')
