from reaction import Reaction
from utils import *
class RecogData():
    def __init__(self):
        self.ocr_list = []
        self.smiles_list = []
        self.ocr_positions = []
        self.smiles_positions = []
    def get(self):
        print(1)
class ChemData():
    def __init__(self):
        self._raw_data = RecogData()
    @property
    def raw_data(self):
        return self._raw_data
    def process_raw_data(self, ):
        self.possible_compounds = ChemUtils.search_compounds_list(self.raw_data.ocr_list)
        self.possible_enzymes = ChemUtils.search_enzymes_list(self.raw_data.ocr_list)
    
    
if __name__ == "__main__":
    a = ChemData()
    a.raw_data.get()