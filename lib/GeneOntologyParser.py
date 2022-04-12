# *_* coding: utf-8 *_*
# @Time     : 11/18/2021 12:54 PM
# @Author   : zong hui
# @object   : some words here


"""
1. we downloaded file go.obo from Gene Ontology website, which is released in 2021-05-01.
2. extract GO information, saved as *../knol/GO/go.info.csv*.
3. extract immune-GO information, saved as *../knol/GO/immunego.info.csv*.
4. GO.obo introduction: https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_4.html 
"""


import pandas as pd
from tqdm import tqdm
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag


class GOPaser:
    """
    - For each GO term, we extract id, name, namespace, level, depth, alt id, 
      ancestors, parents, children, all children, definition
    - We used tools goatool.
    - we obtained all children GO terms of GO:0002376(immune system process), and take it as immune related GO.
      term GO:0002376 (immune system process) is a child term of GO:0008150 (biological process)
    """
    def __init__(self, GOobo):
        self.GOobo = GOobo
        self.godag = GODag(self.GOobo)
        self.go_ids = sorted(list(set([go_term.id for go_term in self.godag.values()])))

    def get_name(self, go_id):
        try:
            return self.godag.get(go_id).name
        except:
            return ''

    def get_namespace(self, go_id):
        try:
            return self.godag.get(go_id).namespace
        except:
            return ''

    def get_level(self, go_id):
        try:
            return self.godag.get(go_id).level
        except:
            return ''

    def get_depth(self, go_id):
        try:
            return self.godag.get(go_id).depth
        except:
            return ''

    def get_alt_ids(self, go_id):
        try:
            return ';'.join(list(self.godag.get(go_id).alt_ids))
        except:
            return ''

    def get_ancestors(self, go_id):
        try:
            gosubdag_r0 = GoSubDag([go_id], self.godag, prt=None)
            return ';'.join(list(gosubdag_r0.rcntobj.go2ancestors[go_id]))
        except:
            pass

    def get_parents(self, go_id):
        try:
            return ';'.join([term.id for term in self.godag.get(go_id).parents])
        except:
            return ''

    def get_children(self, go_id):
        try:
            return ';'.join([term.id for term in self.godag.get(go_id).children])
        except:
            return ''

    def get_all_children(self, go_id):
        try:
            return sorted(list(self.godag.get(go_id).get_all_children()))
        except:
            return ''

    def get_definition(self, go_id):
        return ''

    def get_info(self, go_id):
        info = {'id': go_id, 
                'name': self.get_name(go_id),
                'namespace': self.get_namespace(go_id),
                'level': self.get_level(go_id), 
                'depth': self.get_depth(go_id), 
                'alt_ids': self.get_alt_ids(go_id), 
                'ancestors': self.get_ancestors(go_id), 
                'parents': self.get_parents(go_id),
                'children': self.get_children(go_id),
                'definition': self.get_definition(go_id),
                }
        return info


    def save_info(self, go_ids, save_path=None):
        infos = []
        for go_id in tqdm(go_ids, ncols=80):
            info = self.get_info(go_id)
            infos.append(info)
        df = pd.DataFrame(infos)
        df.to_csv(save_path, index=False)
                

class GOMapping:
    def __init__(self, GOinfo):
        self.GOinfo = pd.read_csv(GOinfo, keep_default_na=False)
        
    def get_id2name(self):
        id2name = dict(zip(self.GOinfo.id, self.GOinfo.name))
        return id2name

    def get_id2namespace(self):
        id2namespace = dict(zip(self.GOinfo.id, self.GOinfo.namespace))
        return id2namespace

    def get_name2id(self):
        name2id = dict(zip(self.GOinfo.name, self.GOinfo.id))
        return name2id

    def get_altid2id(self):
        altid2id = dict()
        for go_id, alt_ids in zip(self.GOinfo.id, self.GOinfo.alt_ids):
            if alt_ids:
                alt_ids_list = alt_ids.split(';')
                for alt_id in alt_ids_list:
                    altid2id[alt_id] = go_id
        return altid2id

    def get_id2level(self):
        id2level = dict(zip(self.GOinfo.id, self.GOinfo.level))
        return id2level

    def get_id2depth(self):
        id2depth = dict(zip(self.GOinfo.id, self.GOinfo.depth))
        return id2depth

    def get_id2ancestors(self):
        id2ancestors = dict(zip(self.GOinfo.id, self.GOinfo.ancestors))
        return id2ancestors

    def get_id2parents(self):
        id2parents = dict(zip(self.GOinfo.id, self.GOinfo.parents))
        return id2parents

    def get_id2children(self):
        id2children = dict(zip(self.GOinfo.id, self.GOinfo.children))
        return id2children

    def get_id2definition(self):
        id2definition = dict(zip(self.GOinfo.id, self.GOinfo.definition))
        return id2definition



if __name__ == "__main__":
    # ====================
    # GOP = GOPaser('../knol/GO/go.obo')
    # get and save information of all GO IDs
    # go_ids = GOP.go_ids
    # print(len(go_ids))
    # GOP.save_info(go_ids, '../knol/GO/go.info.csv')
    # get and save information of all children of GO:0002376 (immune system process)
    # go_id = "GO:0002376"
    # go_id_all_children = GOP.get_all_children('GO:0002376')
    # go_ids = [go_id] + go_id_all_children
    # GOP.save_info(go_ids, '../knol/GO/immunego.info.csv')

    # # ====================
    # GOM = GOMapping('../knol/GO/go.info.csv')
    # name2id = GOM.get_name2id()
    # id2name = GOM.get_id2name()
    # id2namespace = GOM.get_id2namespace()
    # id2level = GOM.get_id2level()
    # id2depth = GOM.get_id2depth()
    # id2ancestors = GOM.get_id2ancestors()
    # id2parents = GOM.get_id2parents()
    # id2children = GOM.get_id2children()
    # id2definition = GOM.get_id2definition()
    # print('id:', name2id.get('reproduction'))
    # print('name:', id2name.get('GO:0000003'))
    # print('namespace:', id2namespace.get('GO:0000003'))
    # print('level:', id2level.get('GO:0000003'))
    # print('depth:', id2depth.get('GO:0000003'))
    # print('ancestors:', id2ancestors.get('GO:0000003'))
    # print('parents:', id2parents.get('GO:0000003'))
    # print('children:', id2children.get('GO:0000003'))
    # print('definition:', id2definition.get('GO:0000003'))
    # altid2id = GOM.get_altid2id()
    # print(altid2id.get('GO:0044152'))
    

    print('done!')
