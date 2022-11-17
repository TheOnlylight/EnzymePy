import pubchempy as pcp
import pickle
import json
from tqdm import tqdm
from rdkit import DataStructs, Chem
class Enzyme():
    def __init__(
        self,
        standard_name = '',
        ec_id =  -1,
):
        self.standard_name = standard_name
    @classmethod
    def dissolve_synonym():
        pass
        
class Compound():
    def __init__(
        self,
        input,
        init_mode
    ):
        try:
            self.pcp_data = pcp.get_compounds(input, init_mode)
            self.smiles = self.pcp_data.isometric_smiles
            self.cid = self.pcp_data.cid
            self.pcp_valid = True
        except Exception:
            self.pcp_valid = False
        
        self.rd_valid = True
        if self.pcp_valid:
            self.rd_data = Chem.MolFromSmiles(self.smiles)
        elif init_mode == 'smiles':
            self.rd_data = Chem.MolFromSmiles(input)
        else:
            self.rd_valid = False
            
        
    def calc_similarity(self, other):
        if self.rd_valid and other.rd_valid:
            ms = [self.rd_data, other.rd_data]
            fps = [Chem.RDKFingerprint(x) for x in ms]
            return DataStructs.FingerprintSimilarity(fps[0],fps[1])
        else:
            return 0
class Reaction():
    def __init__(self):
        self._enzyme = Enzyme()
    @property
    def enzyme(self):
        return self._enzyme