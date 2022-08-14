# *_* coding: utf-8 *_*
# @Time     : 07/01/2021
# @Author   : zong hui
# @object   : some configuration of GOF.


import os


class basicConfig(object):
    # GO_INFO = "../knol/GO/immunego.info.csv"
    GO_INFO = "../knol/GO/go.info.csv"
    DB_PATH = 'BiologicalEntity'
    # TCGA_RNAseq_PATH = '../../rawdata/TCGA_RNAseq/'
    MSigDB_C5_GO_PATH = "../knol/MSigDB/c5.go.v7.2.symbols.gmt"
    NCBI_Entrez_Gene2GO = "../knol/NCBI_EntrezGene/gene2go.gz"
    HumanGeneInformation = "../knol/HumanGeneInformation/HumanGeneInformation.txt"

    
class Config(basicConfig):
    def __init__(self, case, pmid):
        
        self.pmid = pmid
        self.case_name = case
        self.case_path = os.path.split(pmid)[0]

        self.logFile = os.path.join(self.case_path, '{}.log'.format(case))
        
        self.pmid2go = '../case/{}/{}/{}_pubmed2go.csv'.format(case, basicConfig.DB_PATH, case)
        self.pmid2gene = '../case/{}/{}/{}_pubmed2gene.csv'.format(case, basicConfig.DB_PATH, case)

        self.mapping_go2pmid = os.path.join(self.case_path, '{}@mapping_go2pmid.txt'.format(case))
        self.mapping_gene2pmid = os.path.join(self.case_path, '{}@mapping_gene2pmid.txt'.format(case))
        self.mapping_pmid2entities = os.path.join(self.case_path, '{}@mapping_pmid2entities.txt'.format(case))

        self.GOF = os.path.join(self.case_path, '{}@GOF.csv'.format(case))
        self.GGSS = os.path.join(self.case_path, '{}@GGSS.csv'.format(case))
        self.GGSS001 = os.path.join(self.case_path, '{}@GGSS001.csv'.format(case))
        
        # self.GEX =  os.path.join(basicConfig.TCGA_RNAseq_PATH, '{}__gene.normalized_RNAseq__tissueTypeAll.txt'.format(case))
        # self.matirx_GOF = os.path.join(self.case_path, '{}@matirx_GOF.txt'.format(case))
        # self.matrix_GEX = os.path.join(self.case_path, '{}@matrix_GEX.txt'.format(case))
        # self.matrix_GOFxGEX = os.path.join(self.case_path, '{}@matrix_GOFxGEX.txt'.format(case))
        
        # self.GOF_gmt = '../case/{}/GOF_gmt/'.format(case)
        # self.MSigDB_folder = '../case/{}/MSigDB/'.format(case)

