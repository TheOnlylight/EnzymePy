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
    def predict_reactions(self,):
        self.only_enzyme = [ChemUtils.find_reaction(x,) for x in self.possible_enzymes]
        self.only_cid = [ChemUtils.find_reaction(ec = [], cid=[x.cid]) for x in self.possible_compounds if x.pcp_valid]
    
if __name__ == "__main__":
    a = ChemData(
        ['Glycerol', 'NADP', 'GlyDH', 'Glyceraldehyde'],
        ['C(C(CO)O)O']
    )
    a.process_raw_data()
    print(a.possible_enzymes)
    a.predict_reactions()
    print(a.only_enzyme)
    print(a.only_cid)