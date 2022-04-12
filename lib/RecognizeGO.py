# *_* coding: utf-8 *_*
# @Time     : 10/21/2021 2:45 PM
# @Author   : zong hui
# @object   : some words here


import re
import json
import pandas as pd
from tqdm import tqdm


def exact_string_match(terms: list, texts: list) -> list:
    """使用全字符匹配识别文本中的GO术语,
    1.检测词边界
    2.最长字符串匹配原则
    """
    # 将术语按长度从大到小排列
    terms.sort(key=lambda i: len(i), reverse=True)
    terms_lower = [term.lower() for term in terms]
    matched_terms = list()
    # 循环处理所有文本
    for i, sentence in enumerate(tqdm(texts, total=len(texts), ncols=80)):
        sentence_lower = sentence.lower()
        matched_term = list()
        # 循环处理每个术语
        for term,term_lower in zip(terms, terms_lower):
            # 先判断术语在不在文本里
            if term_lower in sentence_lower:
                try:
                    # \\b为词边界，如果匹配到term，则替换掉
                    new_sentence = re.sub('\\b' + term_lower + '\\b', '', sentence_lower)
                except:
                    # 因特殊字符出现re错误，则直接替换
                    new_sentence = sentence_lower.replace(term_lower, '')
                # 新旧句子不一样，说明匹配到了术语
                if new_sentence != sentence_lower:
                    # 更新句子，并添加匹配到的术语
                    sentence_lower = new_sentence
                    matched_term.append(term)
        # 保存每条文本识别到的术语
        matched_terms.append(matched_term)
    return matched_terms

if __name__ == '__main__':
    # ======test======
    terms = ['Immune response', 'signaling', 'membrane', 'regulation of cell aging', 'cell aging', 'aging', 'GO:0016020']
    texts = ['TSPAN6 belongs to the tetraspanin superfamily of multi-pass membrane proteins interacting with multiple immune-related molecules, including immune receptors, integrins, signaling molecules and functioning as important immune response modulators.', 'regulation of cell aging']
    rst = exact_string_match(terms, texts)
    print(rst)
    # =====output=====
    # [['Immune response', 'signaling', 'membrane'], ['regulation of cell aging']]
