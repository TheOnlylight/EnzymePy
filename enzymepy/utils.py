import pubchempy as pcp
import pickle
import json
from tqdm import tqdm
from rdkit import DataStructs, Chem
from nltk.metrics import *
def load_json(path):
    with open(path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    return data

class ChemUtils():
    @staticmethod
    def dynamic_leven_dist(str1, str2):
        result = edit_distance_align(str1, str2, substitution_cost=1)
        return result
    @staticmethod
    def leven_dist(str1, str2):
        return edit_distance(str1, str2, substitution_cost=1)
    @staticmethod
    def search_compounds_list(mols):
        ans = []
        for m in mols:
            ans += pcp.get_compounds(m, 'name')
        return ans
    @staticmethod
    def str_sim(str1, str2):
        tot_len = len(str1)
        dist = ChemUtils.leven_dist(str1, str2)
        good = tot_len - dist
        if(good <= 0):
            good = 0
            return good/tot_len
        else:
            return good/tot_len
    @classmethod
    def load_brenda(
        cls, 
        brenda_datapath = 'data/BrendaIDwithCid_duplicated.pkl', 
        syn_datapath = ''
        ):
        with open(brenda_datapath, 'rb') as handle:
            cls.brenda = pickle.load(handle)
        print("brenda length:", len(cls.brenda))
    @classmethod
    def init_syns(cls):    # load all synonyms of ec
        '''
        syns_list: all synonyms of ec
        reverse_dict: map synonyms to standard names
        '''
        syn_file = 'data/synonyms.json'
        syns = load_json(syn_file)
        syns_list = []
        for key in syns:
            syns[key] += [key]
            syns[key] = list(set(syns[key]))
            syns_list += syns[key]
        print(len(syns_list))

        reverse_dict = {}
        for key in syns:
            for item in syns[key]:
                if item.lower() in reverse_dict:
                    reverse_dict[item.lower()].append(key.lower())
                else:
                    reverse_dict[item.lower()] = [key.lower()]
        print(len(reverse_dict))
        cls.reverse_dict = reverse_dict
        cls.dict = syns
        return syns, syns_list, reverse_dict
    @classmethod
    def get_syns(cls, enzyme):
        return cls.dict[enzyme]
    @classmethod
    def dissolve_enzyme_synonym(cls, name):
        try:
            return cls.reverse_dict[name.lower()]
        except Exception:
            print(f'no enzyme found for {name}')
            return []
    @classmethod
    def find_reactions(cls, duplicate_pairs, ):
        brenda = cls.brenda
        reaction_pairs = []
        for item in duplicate_pairs:
            ec = item[0]
            cid = item[1]
            # print(ec, cid)
            reaction_ids = cls.find_reactions(ec, cid)
            if reaction_ids != []:
                reaction_pairs.append([ec, cid, reaction_ids])
        print("reaction pairs:", len(reaction_pairs))
        return reaction_pairs
    @classmethod
    def find_reaction(cls, ec = [], cid = []):
        """find reaction according to ec and cids

        Args:
            ec (_type_): the enzyme standardized name
            cid (_type_): the cid
        """
        brenda = cls.brenda
        # print(ec)
        reaction_ids = [] # save mapped react id and ent name
        if ec != []:
            for key in brenda:
                if brenda[key]['ec_name'].lower() == ec:
                    set_ex_cid = set(cid)
                    for idx, cand_cids in enumerate(brenda[key]['cids']):
                        set_cand = set(cand_cids)
                        # print(set_cand)
                        if (set_cand & set_ex_cid) or cid == []:
                            cand_ent_name = brenda[key]['cems'][idx] # save the mapped ent name
                            reaction_ids.append([key, cand_ent_name])
        else:
            print(cid)
            set_ex_cid = set(cid)
            for key in brenda:
                for idx, cand_cids in enumerate(brenda[key]['cids']):
                    set_cand = set(cand_cids)
                    # print(set_cand)
                    if (set_cand & set_ex_cid) or cid == []:
                        cand_ent_name = brenda[key]['cems'][idx] # save the mapped ent name
                        reaction_ids.append([key, cand_ent_name])
        return reaction_ids
    @classmethod
    def get_brenda_reaction(cls, id):
        return cls.brenda[id]
ChemUtils.load_brenda()
ChemUtils.init_syns()

if __name__ == "__main__":
    ChemUtils.load_brenda()
    ChemUtils.init_syns()
    a = ChemUtils.dissolve_enzyme_synonym('GlyDH')
    print(a)
    print(ChemUtils.leven_dist('gly', '-glyfdas8'))
    print(ChemUtils.str_sim('gly', '-gly'))