import pubchempy as pcp
import pickle
import json
from tqdm import tqdm
from rdkit import DataStructs, Chem
from utils import ChemUtils
class Enzyme():
    def __init__(
        self,
        standard_name = '',
        name =  '',
    ):
        self.standard_name = standard_name
        self.name = name
        self.synonyms = ChemUtils.get_syns(standard_name) if self.standard_name else []
    def calc_similarity(self, other):
        ans = 0
        for syn in other.synonyms:
            ans = max(ans, ChemUtils.str_sim(syn, self.name))
        return ans
class Compound():
    def __init__(
        self,
        input,
        init_mode
    ):
        self.name = input if init_mode == 'name' else None
        try:
            self.pcp_data = pcp.get_compounds(input.lower(), init_mode)[0]
            self.cid = self.pcp_data.cid
            self.compound = pcp.Compound.from_cid(self.cid)
            self.smiles = self.pcp_data.isomeric_smiles
            self.name = self.pcp_data.iupac_name
            self.pcp_valid = True
        except Exception:
            self.pcp_valid = False
            print(f'no found pcp for {input}',)
        
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
        elif self.name:
            return ChemUtils.str_sim(self.name, other.name)
        else:
            return 0
class Reaction():
    def __init__(self, substrate = None, products = None, enzyme = None, data = {}):
        self._enzyme = enzyme if enzyme else Enzyme()
        self.substrate = substrate
        self.products = products
        if data:
            self._enzyme = Enzyme(standard_name = data['ec_name'])
            self.compounds = [Compound(input = x, init_mode='name') for x in data['cems']]
    @property
    def enzyme(self):
        return self._enzyme
    def similarities(self, compounds, enzymes):
        """this is the place for finding similarities

        Args:
            compounds (list): the possible compounds list
            enzymes (list): the possible enzymes
        """
        tot_cnt = 0
        sim_com = 0
        for com in self.compounds:
            tot_cnt += 1
            for c in compounds:
                sim_com += com.calc_similarity(c)
        sim_com /= tot_cnt
        self.sim_compounds = sim_com