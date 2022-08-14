# *_* coding: utf-8 *_*
# @Time     : 12/9/2020 1:41 PM
# @Author   : zong hui
# @object   : main function of GOF


# load packages
import os
import time
import logging
import datetime
import argparse


# load our modules
import config
from EntityMapping import EntityMapping
from CalculatePvalue import CalculatePvalue
from GeneSimilarity import GeneSimilarity
from GOFxGEX import GOFxGEX


# parsing parameters
parser = argparse.ArgumentParser()
parser.add_argument('--case', '-case', help='a cancer name')
parser.add_argument('--pmid', '-pmid', help='a file, each line is a pmid')
args = parser.parse_args()
case = args.case
pmid = args.pmid


config = config.Config(case, pmid)
logging.basicConfig(level=logging.DEBUG, filename=config.logFile, filemode="w", format="%(message)s")


# step1: Entity Mapping
t1 = time.time()
logging.info("{} GOF, logging:\n{}".format(case, datetime.datetime.now()))
em = EntityMapping(config.pmid, config.pmid2gene, config.pmid2go)
em.save_PMIDEntitesMapping(save_path=config.mapping_pmid2entities)
em.save_GeneMapping(save_path=config.mapping_gene2pmid)
em.save_GoMapping(save_path=config.mapping_go2pmid)
GeneMapping = em.GeneMapping
GoMapping = em.GoMapping
logging.info("\n--------------------\n[step1] Entity Mapping")
logging.info("\t[PMID] number of Input PMID: {}".format(len(em.input_pmid_list)))
logging.info("\t       after screening, leave the PMID, which abstract mentions both Gene (Our HumanGeneInformation) and GO (children of immune system process)")
logging.info("\t       {} PMIDs left, saved at {}".format(len(em.pmid_list), config.mapping_pmid2entities))
logging.info("\t[Gene] {} Genes mentioned in the abstracts, saved at {}".format(len(GeneMapping), config.mapping_gene2pmid))
logging.info("\t[GO]   {} GOs mentioned in the abstracts, saved at {}".format(len(GoMapping), config.mapping_go2pmid))
t2 = time.time()
logging.info("\t[time] used time: {} seconds".format(round(t2-t1, 4)))


# step2: Calculate Enrichment Score:  Pvalue and Adjusted-Pvalue
cp = CalculatePvalue(config.mapping_gene2pmid, config.mapping_go2pmid)
cp.save_GOF(GOF_save_path=config.GOF)
logging.info("\n--------------------\n[step2] Calculate Enrichment Score")
logging.info("\t[PMID] number of Pubmed abstracts mentioned both Gene and GO: {}".format(len(cp.pmid)))
logging.info("\t[GOF]  generating GOF...")
logging.info("\t\t[GOF-Gene] {} Genes in GOF".format(len(cp.GOF_Genes)))
logging.info("\t\t[GOF-GO]   {} GOs in GOF".format(len(cp.GOF_GOs)))
logging.info("\t\t[GOF file] GOF saved at {}".format(config.GOF))
t3 = time.time()
logging.info("\t[time] used time: {} seconds".format(round(t3-t2, 4)))


# step3: Calculate Gene-Gene Similarity Score
gs = GeneSimilarity(config.GOF)
nodes, edges = gs.save_GGSS(save_path=config.GGSS, percentage=1)
logging.info("\t[GGSS]    saved at {}, the node number:[{}], the edge number: [{}]".format(config.GGSS, nodes, edges))
nodes001, edges001 = gs.save_GGSS(save_path=config.GGSS001, percentage=0.01)
logging.info("\n--------------------\n[step3] Calculate Gene-Gene Similarity Score")
logging.info("\t[1% GGSS] The top 1% GGSS saved at {}, the node number:{}, the edge number:{}".format(config.GGSS001, nodes001, edges001))
t4 = time.time()
logging.info("\t[time] used time: {} seconds".format(round(t4-t3, 4)))
logging.info("\t[time] totally use time: {} seconds".format(round(t4-t1, 4)))


# step4: calculate GOFxGEX
# gg = GOFxGEX(config.GOF, config.GEX)
# gg.save_GOFxGEX(save_path = config.matrix_GOFxGEX)
# logging.info("\n--------------------\n[step4] Calculate GOF x GEX Matirx")
# logging.info("\t[GEX] tumor samples:{}, normal samples:{}".format(len(gg.Samples_tumor), len(gg.Samples_normal)))
# logging.info("\t[GOFxGEX matrix: GOxSample] {}x{}, saved at {}".format(len(gg.GOFxGEX_matrix.index), len(gg.GOFxGEX_matrix.columns), config.matrix_GOFxGEX))
# t5 = time.time()
# logging.info("\t[time] used time: {} seconds".format(round(t5-t4, 4)))
# logging.info("\t[time] totally use time: {} seconds".format(round(t5-t1, 4)))

