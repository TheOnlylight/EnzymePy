from .reaction import Reaction,Compound
from .utils import *
import itertools
from tqdm import tqdm
class RecogData():
    def __init__(self,
                ocr_list = [],
                smiles_list = []
                ):
        self.ocr_list = ocr_list
        self.smiles_list = smiles_list
        self.ocr_positions = []
        self.smiles_positions = []
    def get(self):
        print(1)

class ChemData():
    def __init__(self, ocr, smiles):
        self._raw_data = RecogData(ocr, smiles)
        self.correct_data = None
    @property
    def raw_data(self):
        return self._raw_data
    def construct_reaction(self, enzyme_name, substrate, product):
        actual_substrate = []
        actual_product = []
        for j in substrate:
            if j in self.compounds_mapping:
                actual_substrate += [self.compounds_mapping[j]]
            else:
                actual_substrate += j
        for j in product:
            if j in self.compounds_mapping:
                actual_product += [self.compounds_mapping[j]]
            else:
                actual_product += j
        self.add_reaction = Reaction(substrate=actual_substrate, products=actual_product,enzyme=enzyme_name)
        return self.add_reaction
    def process_raw_data(self, strict = False ):
        self.possible_enzymes = []
        self.enzyme_mapping = {}
        self.possible_compounds = []
        self.compounds_mapping = {} # avoid the function
        for e in self.raw_data.ocr_list:
            self.possible_enzymes += ChemUtils.dissolve_enzyme_synonym(e)
            self.possible_enzymes = list(set(self.possible_enzymes))
        for e in self.raw_data.ocr_list:
            # TODO about the input methods and searchings
            x = Compound(input = e, init_mode='name')
            if strict is True and x.pcp_valid is False:
                continue
            self.possible_compounds += [x]
            self.compounds_mapping[e] = x
        for e in self.raw_data.smiles_list:
            x = Compound(input = e, init_mode='smiles')
            self.possible_compounds += [x]
            self.compounds_mapping[e] = x
            
    def predict_pairs(self,):
        self.pairs = []
        valid_cids = [x.cid for x in self.possible_compounds if x.cid]
        self.valid_compounds = [x for x in self.possible_compounds if x.cid]
        
        self.valid_cids = valid_cids
        
        # valid_cids = filter_list(valid_cids, ban_list)
        for x in tqdm(self.possible_enzymes, 'search enzyme'):
            self.pairs+=ChemUtils.find_pairs(x,valid_cids)
    def predict_reactions(self, gross = True, valve = 2, ban_list = []):
        def filter_list(list, ban_list):
            for b in ban_list:
                list = list.remove(b)
            return list
        if gross:
            self.only_enzyme = [ChemUtils.find_reaction(x,) for x in self.possible_enzymes]
            self.only_cid = [ChemUtils.find_reaction(ec = [], cid=[x.cid]) for x in self.possible_compounds if x.pcp_valid]
            self.only_enzyme_reaction = [Reaction(data = ChemUtils.get_brenda_reaction(id[0][0])) for id in self.only_enzyme if id]
            self.only_cid_reaction = [Reaction(data = ChemUtils.get_brenda_reaction(id[0][0])) for id in self.only_cid if id]
        else:
            self.only_enzyme = [ChemUtils.find_reaction(x,) for x in tqdm(self.possible_enzymes, 'search enzyme')]
            valid_cids = [x.cid for x in self.possible_compounds if x.cid]
            self.valid_compounds = [x for x in self.possible_compounds if x.cid]
            
            self.valid_cids = valid_cids
            
            valid_cids = filter_list(valid_cids, ban_list)
            self.valid_reaction = []
            brenda_ids = []
            for e in self.only_enzyme:
                for it in e:
                    brenda_ids.append(it[0])
            brenda_ids = list(set(brenda_ids))
            if valid_cids is not []:
                for id in tqdm(brenda_ids, 'process sub enzyme proc'):
                    # id represent a reaction
                    cur_data = ChemUtils.get_brenda_reaction(id)
                    cur_cids = cur_data['cids']
                    merged_cids = list(itertools.chain(*cur_cids))
                    if len(set(merged_cids)&set(valid_cids)) >= valve:
                        self.valid_reaction += [Reaction(data = cur_data)]
            self.only_enzyme_reaction = [Reaction(data = ChemUtils.get_brenda_reaction(id[0][0])) for id in self.only_enzyme if id]
            self.only_cid_reaction = []

    def show_sim(self):
        for j in self.only_enzyme_reaction:
            j.similarities(compounds=self.possible_compounds, enzymes=self.possible_enzymes)
            print(j.sim_compounds)
        for j in self.only_cid_reaction:
            j.similarities(compounds=self.possible_compounds, enzymes=self.possible_enzymes)
            print(j.sim_compounds)
    def calc_sim(self, calc_method = 'name'):

        cur_sim = 0
        for idx, re in enumerate(self.only_cid_reaction):
            self.only_cid_reaction[idx].similarities(
                compounds=self.possible_compounds,
                enzymes=self.possible_enzymes
            )
            if cur_sim < self.only_cid_reaction[idx].sim_compounds:
                self.bet_ans = self.only_cid_reaction[idx]
                cur_sim = self.only_cid_reaction[idx].sim_compounds
        for idx, re in enumerate(tqdm(self.only_enzyme_reaction, 'calc sims')):
            self.only_enzyme_reaction[idx].similarities(
                compounds=self.possible_compounds,
                enzymes=self.possible_enzymes
            )
            if cur_sim < self.only_enzyme_reaction[idx].sim_compounds:
                self.bet_ans = self.only_enzyme_reaction[idx]
                cur_sim = self.only_enzyme_reaction[idx].sim_compounds
        # self.bet_ans.pprint()
    
if __name__ == "__main__":
    a = ChemData(
        ['Glycerol', 'NADP', 'GlyDH', 'Glycerone'],
        ['C(C(CO)O)O']
    )
    a.process_raw_data()
    # print(a.possible_enzymes)
    a.predict_reactions()
    print(a.only_enzyme[0][0])
    print(ChemUtils.get_brenda_reaction(265)['cems'])
    print(a.only_cid_reaction)
    a.calc_sim()