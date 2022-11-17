import pubchempy as pcp
import pickle
import json
from tqdm import tqdm
from rdkit import DataStructs, Chem

def load_json(path):
    with open(path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    return data

class ChemUtils():
    @staticmethod
    def search_compounds_list(mols):
        ans = []
        for m in mols:
            ans += pcp.get_compounds(m, 'name')
        return ans
    @classmethod
    def load_brenda(
        cls, 
        brenda_datapath = '/data2/private/zhouhantao/MultiModal/kara/data/BrendaIDwithCid_duplicated.pkl', 
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
        syn_file = '/data2/private/zhouhantao/data/synonyms.json'
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
                if item in reverse_dict:
                    reverse_dict[item.lower()].append(key.lower())
                else:
                    reverse_dict[item.lower()] = [key.lower()]
        print(len(reverse_dict))
        cls.reverse_dict = reverse_dict
        return syns, syns_list, reverse_dict
    @classmethod
    def dissolve_enzyme_synonym(cls, name):
        return cls.reverse_dict[name.lower()]
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
    def find_reaction(cls, ec, cid = []):
        """find reaction according to ec and cids

        Args:
            ec (_type_): the enzyme standardized name
            cid (_type_): the cid
        """
        brenda = cls.brenda
        reaction_ids = [] # save mapped react id and ent name
        for key in brenda:
            if brenda[key]['ec_name'] == ec:
                set_ex_cid = set(cid)
                for idx, cand_cids in enumerate(brenda[key]['cids']):
                    set_cand = set(cand_cids)
                    if (set_cand & set_ex_cid) or cid == []:
                        cand_ent_name = brenda[key]['cems'][idx] # save the mapped ent name
                        reaction_ids.append([key, cand_ent_name])


if __name__ == "__main__":
    ChemUtils.load_brenda()
    ChemUtils.init_syns()
    a = ChemUtils.dissolve_enzyme_synonym('GlyDH')
    print(a)