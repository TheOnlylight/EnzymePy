import pubchempy as pcp
import pickle
import json
from tqdm import tqdm
from rdkit import DataStructs, Chem
from .utils import ChemUtils
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
        init_mode,
        init_strict = False,
    ):
        self.name = input if init_mode == 'name' else None
        self.cid = 0
        try:
            try:
                cid = ChemUtils.local_pcp[input]
                if cid == []:
                    self.cid = 0
                else:
                    self.cid = cid[0]
                    # simply take the first cid matching from pubchem
            except:
                self.pcp_data = pcp.get_compounds(input.lower() if init_mode == 'name' else input, init_mode)[0]
                print('online pcp')
                self.cid = self.pcp_data.cid
            # self.compound = pcp.Compound.from_cid(self.cid)
            # self.smiles = self.pcp_data.isomeric_smiles
            # self.name = self.pcp_data.iupac_name
            self.pcp_valid = False
        except Exception:
            self.pcp_valid = False
            # print(f'no found pcp for {input}',)
        
        self.rd_valid = True
        if self.pcp_valid:
            self.rd_data = Chem.MolFromSmiles(self.smiles)
            self.smiles = Chem.MolToSmiles(self.rd_data,kekuleSmiles=True)
        elif init_mode == 'smiles':
            self.rd_data = Chem.MolFromSmiles(input)
            self.smiles = Chem.MolToSmiles(self.rd_data,kekuleSmiles=True)
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
    def show_stucture(self):
        self.image = None
        if self.rd_valid:
            self.image = Chem.Draw.MolToImage(self.rd_data)
            return self.image
class Reaction():
    def __init__(self, substrate = None, products = None, enzyme = None, data = {}):
        self._enzyme = enzyme if enzyme else Enzyme()
        self.substrate = substrate
        self.products = products
        

        if data:
            self.cems = data['cems']
            self.cids = data['cids']
            self._enzyme = Enzyme(standard_name = data['ec_name'])
            # self.compounds = [Compound(input = x[0], init_mode='cid') for x in data['cids'] if x]
            self.compounds = []
            self.data = data
            pivot = data['react'].find('=')
            substrate = []
            products = []
            for j in data['cems']:
                if data['react'].find(j) < pivot:
                    substrate += [j]
                else:
                    products += [j]
            self.substrate = substrate
            self.products = products

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
                sim_com = max(sim_com, com.calc_similarity(c))
        sim_com /= tot_cnt if tot_cnt > 0 else 1
        self.sim_compounds = sim_com
    def pprint(self,):
        print(self.enzyme, self.compounds)
        print(self.cems)
    def get_images(self):
        self.images = []
        for j in self.compounds:
            self.images += [j.show_stucture()]
    def get_dict(self):
        ans = {}
        ans['ec_name'] = self.enzyme
        reaction = ""
        cems = []
        for j in self.substrate:
            reaction += j.name if type(j) is not str else j
            cems += [j.name if type(j) is not str else j]
        reaction += " = "
        for j in self.products:
            reaction += j.name if type(j) is not str else j
            cems += [j.name if type(j) is not str else j]
        ans['react'] = reaction
        ans['cems'] = cems
        cids = []
        for j in cems:
            try:
                cids += [pcp.get_cids(j)]
            # trunk-ignore(flake8/E722)
            except:
                cids += [[]]
        ans['cids'] = cids # TODO
        return ans