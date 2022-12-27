import pubchempy as pcp
import pickle
from rdkit import DataStructs, Chem
from nltk.metrics import *
import pkgutil
import json

data_syn = pkgutil.get_data(__name__, "data/syn.pkl")
data_brenda = pkgutil.get_data(__name__, "data/BrendaIDwithCid_duplicated.pkl")
class ChemUtils():
    @classmethod
    def load_local_cid(cls, file_path):
        with open(file_path, 'r') as f:
            cls.local_pcp = json.load(f)
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
    def add_entry_brenda(cls, entry):
        new_id = list(cls.brenda.keys())[-1] + 1
        cls.brenda[new_id] = entry
        cls.dict[entry['ec_name'].lower()] = [entry['ec_name'].lower()]
        cls.reverse_dict[entry['ec_name'].lower()] = [entry['ec_name'].lower()]
    @classmethod
    def load_data(cls):
        syns = pickle.loads(data_syn)
        cls.brenda = pickle.loads(data_brenda)
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
        return cls.dict[enzyme.lower()]
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
                if brenda[key]['ec_name'].lower() == ec.lower():
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
ChemUtils.load_data()

if __name__ == "__main__":
    ChemUtils.load_brenda()
    ChemUtils.init_syns()
    a = ChemUtils.dissolve_enzyme_synonym('GlyDH')
    print(a)
    print(ChemUtils.leven_dist('gly', '-glyfdas8'))
    print(ChemUtils.str_sim('gly', '-gly'))