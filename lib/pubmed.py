# *_* coding: utf-8 *_*
# @Time     : 10/11/2021 12:59 PM
# @Author   : zong hui
# @object   : some useful functions related to PubMed.gov database

"""
功能主要包括：
    term2pmid:          输入查询词，访问PubMed数据库，获得PMID列表。
    parsePMCIDsFile:    解析PubMed数据库的ID映射文件，获得PMID和PMCID的映射。
    pmid2pmcid:         输入pmid列表，获取对应的pmcid。
    downloadMedline:    输入pmid列表，下载Medline格式数据。
"""


import os
import csv
import json
import time
import requests
import pandas as pd
from Bio import Medline, Entrez
Entrez.email = ""


def term2pmid(terms: list, save_path=None) -> list:
    """输入查询词，访问PubMed数据库，获得PMID列表."""
    # 构建查询语句
    query_terms = '(' + ') OR ('.join(terms) + ')'
    print('[query]: {}'.format(query_terms))
    # 开始检索
    handle0 = Entrez.esearch(db='pubmed', term=query_terms, RetMax=300000000)
    record = Entrez.read(handle0)
    pmids, count = record['IdList'], record['Count']
    # 按升序进行排序
    pmids.sort(key=int)
    # 保存结果
    if save_path:
        with open(save_path, 'w') as f:
            f.write('\n'.join(pmids) + '\n')
        print('[pmid]: {}, [save path]: {}'.format(len(pmids), save_path))
    # 休眠1秒，避免持续访问导致连接中断。
    time.sleep(1)
    return pmids


def parsePMCIDsFile(PMCIDS_file = "../knol/PubMed/PMC-ids.csv"):
    """
    The PMC-ids.csv.gz file, available through the FTP service,
    maps an article’s standard IDs to each other and to other article metadata elements.
    PMC-ids.csv.gz is a comma-delimited file with the following fields:
    Journal Title, ISSN, ..., PMCID, PubMed ID (if available), ..., Release Date (Mmm DD YYYY or live)
    link: https://www.ncbi.nlm.nih.gov/pmc/pmctopmid/
    """
    PMCIDtoPMID = dict()
    PMIDtoPMCID = dict()
    with open(PMCIDS_file, 'r') as inf:
        rows = list(csv.reader(inf))
        for row in rows[1:]:
            PMCIDtoPMID[row[8]] = row[9]
            if row[9]:
                PMIDtoPMCID[row[9]] = row[8]
    print('{} pmcid mapped to {} pmid in file PMC-ids.csv.gz'.format(len(PMCIDtoPMID), len(PMIDtoPMCID)))
    return PMCIDtoPMID, PMIDtoPMCID


def pmid2pmcid(pmids: list, PMIDtoPMCID: dict, save_path=None):
    """输入PMID，获得相应的PMCID"""
    pmcids = [PMIDtoPMCID.get(p, '') for p in pmids]
    pmcids_count = len(pmcids) - pmcids.count('')
    if save_path:
        with open(save_path, 'w') as f:
            f.write('{}\t{}\n'.format('PMID', 'PMCID'))
            for pmid, pmcid in zip(pmids, pmcids):
                f.write('{}\t{}\n'.format(pmid, pmcid))
        print('[pmid]: {}, [pmcid]: {}, [save path]: {}'.format(len(pmids), pmcids_count, save_path))
    return pmcids


def downloadMedline(case, pmids:list=None, save_path=None):
    t1 = time.time()
    # 如果存储路径不存在，则创建
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    print("[case]: {}, [pmid]: {}, [save path]:{}".format(case, len(pmids), save_path))
    # 分批次下载，每次下载10000条
    count = len(pmids)
    batch_size = 10000
    iterations = [[i * batch_size, min((i + 1) * batch_size, count)] for i in range((count-1) // batch_size + 1)]
    # 开始分批次下载
    for (start, end) in iterations:
        # 获取数据，Medline格式
        handle1 = Entrez.efetch(db='pubmed', id=pmids[start:end], rettype='medline', retmode='text')
        record_medline = Medline.parse(handle1)
        # 保存数据，json格式
        cur_save_path = os.path.join(save_path, "{}({}-{}).json".format(case, start+1, end))
        with open(cur_save_path, "w") as f:
            f.write(json.dumps(list(record_medline), ensure_ascii=False, indent=4))
        print('\t[Downloading]: {}-{}, [saved]: {}'.format(start+1, end, cur_save_path))
        # 休眠1秒，避免持续访问导致连接中断。
        time.sleep(1)
    t2 = time.time()
    print('\t[used time]: {} seconds.'.format(round(t2-t1, 4)))
    return "downloaded!"


if __name__ == "__main__":
    #==================================================
    # PMCIDtoPMID, PMIDtoPMCID = parsePMCIDsFile()
    # 输入9种virus的检索词，获取并保存PMID，获取PMID和PMCID的映射
    # df = pd.read_json('../knol/virus/virus_synonym.json')
    # case_terms = df.set_index('case').to_dict()['synonym']
    # for case, terms in case_terms.items():
    #     print('[case]: ', case)
    #     case_folder = '../case/virus/{}'.format(case)
    #     if not os.path.exists(case_folder): os.mkdir(case_folder)
    #     pmids = term2pmid(terms, save_path=os.path.join(case_folder, '{}.pmid.txt'.format(case)))
    #     pmcids = pmid2pmcid(pmids, PMIDtoPMCID, save_path=os.path.join(case_folder, '{}.pmcid.txt'.format(case)))

    # 输入33种TCGA的检索词，获取并保存PMID。
    # df = pd.read_json('../knol/TCGA/TCGA_synonym.json')
    # case_terms = df.set_index('case').to_dict()['synonym']
    # for case, terms in case_terms.items():
    #     print('[case]: ', case)
    #     case_folder = "../case/TCGA/{}".format(case)
    #     if not os.path.exists(case_folder): os.mkdir(case_folder)
    #     pmids = term2pmid(terms=terms, save_path=os.path.join(case_folder, "{}.pmid.txt".format(case)))
    #     pmcids = pmid2pmcid(pmids, PMIDtoPMCID, save_path=os.path.join(case_folder, '{}.pmcid.txt'.format(case)))

    # 输入6个GIDB的检索词，获取并保存PMID。
    # df = pd.read_json('../knol/GIDB/GIDB_synonym.json')
    # case_terms = df.set_index('case').to_dict()['synonym']
    # for case, terms in case_terms.items():
    #     print('[case]: ', case)
    #     case_folder = "../case/GIDB/{}".format(case)
    #     if not os.path.exists(case_folder): os.mkdir(case_folder)
    #     pmids = term2pmid(terms, save_path=os.path.join(case_folder, "{}.pmid.txt".format(case)))
    #     pmcids = pmid2pmcid(pmids, PMIDtoPMCID, save_path=os.path.join(case_folder, '{}.pmcid.txt'.format(case)))


    # ==================================================
    # 输入PMID，从PubMed下载Medline格式数据。
    # cases = ['AAV2', 'COVID19', 'EBV', 'HBV', 'HIV', 'HPV', 'HTLV1', 'MCV', 'XMRV']
    # cases = ['EBV']
    # for case in cases:
    #     pmid_file = '../case/virus/{}/{}.pmid.txt'.format(case, case)
    #     with open(pmid_file, 'r') as f:
    #         pmids = [line.strip() for line in f]
    #     save_path = '../data/PubMed/{}'.format(case)
    #     downloadMedline(case=case, pmids=pmids, save_path=save_path)

    # ==================================================
    print('done!')
