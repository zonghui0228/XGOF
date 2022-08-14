# *_* coding: utf-8 *_*
# @Time     : 10/15/2021 11:16 AM
# @Author   : zong hui
# @object   : some useful functions related to PubTator.gov database


"""
功能主要包括：
    call_pubtator_api:    调用Pubtator的API下载注释数据。格式为biocjson，输入idList: PMID(摘要)，PMCID(全文)。
    biocjson_to_list:      提取PubTator的biocjson格式里的注释信息，转换为为list，添加原句信息。
    download_pubtator:   输入pmids和pmcids，下载Pubtator的标注。
    split_anno_sent:      仅保留基因注释，并将注释和原句分开保存，节省存储空间。
"""

import os
import sys
import json
import getopt
import requests
from tqdm import tqdm
from collections import defaultdict
from nltk.tokenize.punkt import PunktSentenceTokenizer

# 代理
PROXIES = {'http': None, 'https': None}
# PROXIES = {'http': 'http://165.225.96.34:10015', 'https': 'https://165.225.96.34:10015'}


def call_pubtator_api(idList, input_type, data_type='biocjson') -> list:
    # URL
    url = 'https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/'
    # 实体类型
    concepts = ['gene', 'disease', 'chemical', 'species', 'mutation', 'cellline']
    dj = {'concepts': concepts}
    if input_type == 'pmid':
        dj = {'pmids': idList}
    if input_type == 'pmcid':
        dj = {'pmcids': idList}
    res = requests.post(url=url + data_type, json=dj, proxies=PROXIES)
    annotations = []
    if res.status_code != 200:
        print("\t[Error]: HTTP code " + str(res.status_code))
    else:
        annotations = res.text.strip().split('\n')
    return annotations


def biocjson_to_list(biocjson: dict) -> list:
    """解析biocjson格式的数据，获得其中的标注"""
    # 获得pmid, pmcid, passages
    pmid = biocjson.get('pmid')
    pmcid = biocjson.get('pmcid')
    passages = biocjson.get('passages')
    # 开始解析biocjson格式注释
    results = []
    for passage in passages:
        # 获取当前被注释的文本段
        curr_text = passage.get('text').strip()
        # 将注释文本段分句，包括起始位置、终止位置
        sentences = [[start, end, curr_text[start:end]] for start, end in
                     PunktSentenceTokenizer().span_tokenize(curr_text)]
        # 当前文本段在全文里的起始位置
        text_offset = passage.get('offset')
        # 当前文本段在全文里的类型
        section = passage.get('infons').get('section')
        section_type = passage.get('infons').get('section_type')
        # 当前文本段的注释信息
        annotations = passage.get('annotations')
        for annotation in annotations:
            ann_text = annotation['text']  # 注释字段
            ann_type = annotation['infons']['type']  # 注释字段的类型
            ann_iden = annotation['infons']['identifier']  # 注释字段的标准化结果
            ann_length = len(annotation['text'].strip())  # 注释字段长度
            # ann_length = annotation['locations'][0]['length'] # 注释字段长度
            ann_full_text_start = annotation['locations'][0]['offset']  # 注释字段的全文起始位置
            ann_full_text_end = annotation['locations'][0]['offset'] + ann_length  # 注释字段的全文终止位置
            ann_curr_text_start = ann_full_text_start - text_offset  # 注释字段在当前文本段的起止位置
            ann_curr_text_end = ann_full_text_end - text_offset  # 注释字段在当前文本段的终止位置
            # 如果注释字段的起止位置在当前文本段外，则跳过
            if ann_curr_text_start < 0 or ann_curr_text_end > len(curr_text):
                continue
            # 注释字段的所属句子在当前文本段中起始位置、终止位置
            anchor_sentence_start = [start for start, end, sentence in sentences if start <= ann_curr_text_start][-1]
            anchor_sentence_end = [end for start, end, sentence in sentences if end >= ann_curr_text_end][0]
            # 所属句子
            anchor_sentence = curr_text[anchor_sentence_start:anchor_sentence_end]
            anchor_sentence = anchor_sentence.replace('\r\n', '  ').replace('\n', ' ')
            # 注释字段在所属句子中的起始位置
            anchor_ann_start = ann_curr_text_start - anchor_sentence_start
            # 注释字段在所属句子中的终止位置
            anchor_ann_end = ann_curr_text_end - anchor_sentence_start
            anchor_ann_text = anchor_sentence[anchor_ann_start: anchor_ann_end]
            results.append([pmid, ann_full_text_start, ann_full_text_end, ann_text, ann_type, ann_iden, pmcid, section,
                            section_type, anchor_ann_start, anchor_ann_end, anchor_ann_text, anchor_sentence])
    # 按照在全文的起始位置进行排序
    results = sorted(results, key=lambda x: x[2])
    return results


def download_pubtator(case, idList, input_type, save_folder, batch_save=False):
    """输入pmid或者pmcids，下载、解析并保存来自pubtator的注释"""
    # 分批次下载，每次下载500条
    batch_size = 500
    id_count = len(idList)
    id_iterations = [[i * batch_size, min((i + 1) * batch_size, id_count)] for i in
                     range((id_count - 1) // batch_size + 1)]
    # 每个batch下载的biocjson都保存，主要为避免长时间下载自动中断
    batch_save_boicjson=True
    pubtator_biocjson = []
    for (start, end) in id_iterations:
        if batch_save_boicjson:
            save_path = os.path.join(save_folder, f'{case}.{input_type}.pubtator.biocjson({start+1}-{end}).txt')
            if os.path.exists(save_path):
                print(f'\t[downloaded]: {start + 1} - {end}')
                with open(save_path, 'r', encoding='utf-8') as f:
                    annotations = [line.strip('\n') for line in f]
            else:
                print(f'\t[downloading]: {start + 1} - {end}')
                annotations = call_pubtator_api(idList[start:end], input_type)
                with open(save_path, 'w', encoding='utf-8') as f:
                    f.write('\n'.join(annotations))
                print(f'\t[pubtator biocjson save path]:', save_path)
        else:
            print(f'\t[downloading]: {start + 1} - {end}')
            annotations = call_pubtator_api(idList[start:end], input_type)
        pubtator_biocjson.extend(annotations)
    # 保存biocjson
    save_boicjson = False
    if save_boicjson:
        save_path = os.path.join(save_folder, f'{case}.{input_type}.pubtator.biocjson.txt')
        with open(save_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(pubtator_biocjson))
        print('\t[pubtator biocjson save path]: ', save_path)
    # biocjson转换为bioclist
    pubtator_bioclist = []
    for ann_biocjson in tqdm(pubtator_biocjson):
        ann_bioclist = biocjson_to_list(json.loads(ann_biocjson))
        pubtator_bioclist.extend(ann_bioclist)
    # 保存bioclist
    save_bioclist = False
    if save_bioclist:
        save_path = os.path.join(save_folder, f'{case}.{input_type}.pubtator.bioclist.txt')
        with open(save_path, 'w', encoding='utf-8') as f:
            for ann_sent in pubtator_bioclist:
                f.write('\t'.join([str(a) for a in ann_sent]) + '\n')
        print('\t[pubtator bioclist save path]: {}', save_path)
    return pubtator_bioclist


def split_anno_sent(case, pubtator_bioclist, save_folder):
    """将标注和原句分开保存"""
    # 仅保留Gene注释
    anno_sent = [a for a in pubtator_bioclist if a[4] == 'Gene']
    # 获得所有原句
    pmid2sent = defaultdict(list)
    for l in anno_sent:
        pmid2sent[l[0]].append(l[-1])
    # 每条原句分配一个ID
    pmidsent2id = {(p, s): str(p) + '_s' + str(i) for p, ss in pmid2sent.items() for i, s in
                   enumerate(sorted(list(set(ss))))}
    # 获得所有标注
    annos = []
    for l in anno_sent:
        # 最后一列,由原句改为原句ID
        annos.append(l[:-1] + [pmidsent2id.get((l[0], l[-1]))])
    # 保存结果
    with open(os.path.join(save_folder, f'{case}.pubtator.bioclist.sent.txt'), 'w') as f:
        col_names = ['pmid', 'sentence_id', 'sentence']
        f.write('\t'.join(col_names) + '\n')
        for k, v in pmidsent2id.items():
            f.write(f'{k[0]}\t{v}\t{k[1]}\n')
    with open(os.path.join(save_folder, f'{case}.pubtator.bioclist.anno.txt'), 'w') as f:
        col_names = ['pmid', 'fulltext_start', 'fulltext_end', 'text', 'type', 'identifier', 'pmcid', 'section',
                     'section_type', 'sentence_start', 'sentence_end', 'sentence_text', 'sentence_id']
        f.write('\t'.join(col_names) + '\n')
        for anno in annos:
            anno = [str(a) for a in anno]
            f.write('\t'.join(anno) + '\n')
    return "done!"


def get_pmid_pmcid(infile):
    total, pmids, pmcids = [], [], []
    with open(infile, 'r') as f:
        for line in f:
            l = line.strip('\n').split('\t')
            total.append(l[0])
            if l[1]:
                pmcids.append(l[1])
            else:
                pmids.append(l[0])
    return total, pmids, pmcids


def save_sentenceid(infile, outfile):
    with open(infile, 'r') as inf:
        inf.readline()
        sentence_id = [line.strip().split('\t')[1] for line in inf]
    with open(outfile, 'w') as outf:
        outf.write('\n'.join(sentence_id) + '\n')


def main(argv):
    # parse parameters
    case = ''
    try:
        opts, args = getopt.getopt(argv, "hc:", ["case="])
    except getopt.GetoptError:
        print('python pubtator.py -c <case>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python pubtator.py -c <case>')
            sys.exit()
        elif opt in ("-c", "--case"):
            case = arg
    
    # obtain annotation
    pmcid_file = f'../case/{case}/{case}.pmcid.txt'
    total, pmids, pmcids = get_pmid_pmcid(pmcid_file)
    print('[case]:', case, '[total]:', len(total), '[pmid-abstract]:', len(pmids), '[pmcid-fulltetx]:', len(pmcids))
    save_folder = f'../case/{case}/PubTator/'
    if not os.path.exists(save_folder):
        os.mkdir(save_folder)
    pmid_pubtator_bioclist = download_pubtator(case=case, idList=pmids, input_type='pmid', save_folder=save_folder)
    pmcid_pubtator_bioclist = download_pubtator(case=case, idList=pmcids, input_type='pmcid', save_folder=save_folder)
    split_anno_sent(case, pmid_pubtator_bioclist+pmcid_pubtator_bioclist, save_folder)
    save_sentenceid(f'../case/{case}/PubTator/{case}.pubtator.bioclist.sent.txt', f'../case/{case}/{case}.sentid.txt')

if __name__ == "__main__":
    main(sys.argv[1:])

