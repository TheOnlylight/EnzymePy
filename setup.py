#!/usr/bin/env python

from distutils.core import setup

setup(name='enzymepy',
      version='1.0',
      description='Python Tool For Enzyme Processing',
      author='Hantao Zhou',
      author_email='teddy_cou@outlook.com',
      packages=['enzymepy'],
      package_dir = {'enzymepy': 'src/enzymepy'},
      package_data = {'enzymepy':['data/BrendaIDwithCid_duplicated.pkl', 'data/syn.pkl']},
      requires=['nltk', 'tqdm', 'rdkit', 'pubchempy'],
      include_package_data = True,
     )