import os
import json
import pickle
# def load_brenda(
#     brenda_datapath = './BrendaIDwithCid_duplicated.pkl', 
#     syn_datapath = ''
#     ):
#     print(os.getcwd())
#     with open(brenda_datapath, 'rb') as handle:
#         brenda = pickle.load(handle)
#     print("brenda length:", len(brenda))
#     return brenda
# brenda = load_brenda()


# def init_syns():    # load all synonyms of ec
#     '''
#     syns_list: all synonyms of ec
#     reverse_dict: map synonyms to standard names
#     '''
#     syn_file = 'synonyms.json'
#     syns = load_json(syn_file)
#     syns_list = []
#     for key in syns:
#         syns[key] += [key]
#         syns[key] = list(set(syns[key]))
#         syns_list += syns[key]
#     print(len(syns_list))

#     reverse_dict = {}
#     for key in syns:
#         for item in syns[key]:
#             if item.lower() in reverse_dict:
#                 reverse_dict[item.lower()].append(key.lower())
#             else:
#                 reverse_dict[item.lower()] = [key.lower()]
#     print(len(reverse_dict))
#     reverse_dict = reverse_dict
#     return syns, syns_list, reverse_dict
import importlib.resources
importlib.resources.files('syn.pkl', 'BrendaIDwithCid_duplicated.pkl')