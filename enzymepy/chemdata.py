from reaction import Reaction,Compound
from utils import *
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
    def process_raw_data(self, ):
        self.possible_enzymes = []
        self.possible_compounds = []
        for e in self.raw_data.ocr_list:
            self.possible_enzymes += ChemUtils.dissolve_enzyme_synonym(e)
        for e in self.raw_data.ocr_list:
            self.possible_compounds += [Compound(input = e, init_mode='name')]
        for e in self.raw_data.smiles_list:
            self.possible_compounds += [Compound(input = e, init_mode='smiles')]
    def predict_reactions(self,):
        self.only_enzyme = [ChemUtils.find_reaction(x,) for x in self.possible_enzymes]
        self.only_cid = [ChemUtils.find_reaction(ec = [], cid=[x.cid]) for x in self.possible_compounds if x.pcp_valid]
        self.only_enzyme_reaction = [Reaction(data = ChemUtils.get_brenda_reaction(id[0][0])) for id in self.only_enzyme if id]
        self.only_cid_reaction = [Reaction(data = ChemUtils.get_brenda_reaction(id[0][0])) for id in self.only_cid if id]
    def show_sim(self):
        for j in self.only_enzyme_reaction:
            j.similarities(compounds=self.possible_compounds, enzymes=self.possible_enzymes)
            print(j.sim_compounds)
        for j in self.only_cid_reaction:
            j.similarities(compounds=self.possible_compounds, enzymes=self.possible_enzymes)
            print(j.sim_compounds)
    
if __name__ == "__main__":
    a = ChemData(
        ['Glycerol', 'NADP', 'GlyDH', 'Glyceraldehyde'],
        ['C(C(CO)O)O']
    )
    a.process_raw_data()
    # print(a.possible_enzymes)
    a.predict_reactions()
    print(a.only_enzyme[0][0])
    print(ChemUtils.get_brenda_reaction(265)['cems'])
    print(a.only_cid_reaction)
    a.show_sim()