# *_* coding: utf-8 *_*
# @Time     : 11/11/2020 12:59 PM
# @Author   : zong hui
# @object   : given a list of pmid, get the mapping of pmid to gene


import os
import csv
import json
import gzip
import pandas as pd
from tqdm import tqdm
from collections import defaultdict


GOinfos_file = '../knol/GO/go.info.tsv'
GOinfos_df = pd.read_csv(GOinfos_file, sep='\t')
go_idname2id = dict()
for go_id, go_name in zip(GOinfos_df['id'], GOinfos_df['name']):
    go_idname2id[go_id] = go_id
    go_idname2id[go_name] = go_id
    

def getGOWithExactStringMatch(case, pmid_path, pmid2go_path):
    # 获取pmid列表
    with open(pmid_path, 'r') as f:
        pmid = [line.strip() for line in f]
        
    # 基于全字符匹配获得摘要中出现的GO术语
    pmid2go_exist = False # 用来控制是否重新获得GO

    if pmid2go_exist == True:
        print("GO has been recognized")
        gos = list()
        with open(pmid2go_path, 'r') as f:
            rows = list(csv.reader(f))
            for row in rows[1:]:
                gos.extend(row[1].split(';'))
        gos = list(set(gos))
        gos.remove('')
    
    if pmid2go_exist == False:
        print("start recognize GO")
        # 获取摘要中出现的GO术语
        pmid2go = defaultdict(list)
        # 下载的摘要路径
        pubmed_path = "../../rawdata/PubMed/{}/".format(case)
        abstracts_files = [file for file in os.listdir(pubmed_path) if file.endswith('.json')]
        for abstracts_file in abstracts_files:
            # 读取每个摘要,识别其中的GO术语
            with open(os.path.join(pubmed_path, abstracts_file), 'r') as f:
                data = json.load(f)
                for d in tqdm(data):
                    pmid_, abstract = d.get('PMID', ''), d.get('AB', '')
                    if pmid_ and abstract:
                        go_ids = list()
                        # 循环所有GO术语，在摘要中全字符匹配。
                        for i in go_idname2id.keys():
                            if i.lower() in abstract.lower():
                                go_ids.append(go_idname2id.get(i))
                        pmid2go[pmid_] = go_ids

        # 存储摘要中出现的GO术语
        gos = list()
        with open(pmid2go_path, 'w', newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["PubMed ID", "GO Terms"])
            for p in pmid:
                writer.writerow([p, ";".join(pmid2go.get(p,''))])
                gos.extend(pmid2go.get(p, []))
    print(len(pmid), len(gos))
    
if __name__ == "__main__":
    case = 'HIV'
    pmid_path = "../case/case-virus/HIV/HIV.pmid.txt"
    pmid2go_path = "../db/GO/HIV_pubmed2go.csv"
    getGOWithExactStringMatch(case, pmid_path, pmid2go_path)

