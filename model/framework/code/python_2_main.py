#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tues Jul 26 2022

@author: Samuel Volk
"""
from sklearn.externals import joblib
import numpy as np
import pandas as pd
from rdkit.Chem import DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit import Chem
import sys
import os

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
file_path = os.path.dirname(os.path.abspath(__file__))

# hide gpus
os.environ["CUDA_VISIBLE_DEVICES"]="-1"

# point to the admetlab paths
path_root = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "..")
sys.path.append(os.path.abspath(path_root))
sys.path.append(os.path.join(path_root, "checkpoints", "classification_model"))

# configure relevant resource paths
checkpoints_dir = os.path.abspath(os.path.join(file_path, "..", "..", "checkpoints"))
framework_dir = os.path.abspath(os.path.join(file_path, ".."))

# read smiles input and model metadata csv files
smiles_df = pd.read_csv(input_file)
smiles_list = smiles_df[smiles_df.columns[0]].tolist()
model_meta = pd.read_csv(os.path.join(framework_dir, 'model_metadata.csv'))

# function for creating ECFP descriptors from SMILES strings
def smiles_to_ECFP(smiles_list, diameter = 4, numBits = 1024):
    df = pd.DataFrame(columns=[str(i) for i in range(numBits)])
    for i in smiles_list:
        arr = np.zeros((0, ), dtype=np.int8)
        mol = Chem.MolFromSmiles(i)
        ECFP_val = AllChem.GetMorganFingerprintAsBitVect(mol, diameter/2, numBits) 
        DataStructs.ConvertToNumpyArray(ECFP_val, arr)
        df.loc[len(df)] = arr 
    return df.to_numpy()

# function for creating MACCS descriptors from SMILES strings
def smiles_to_MACCS(smiles_list):
    df = pd.DataFrame(columns=[str(i) for i in range(167)])
    for i in smiles_list:
        arr = np.zeros((0, ), dtype=np.int8)
        mol = Chem.MolFromSmiles(i)
        MACCS_val = rdMolDescriptors.GetMACCSKeysFingerprint(mol) 
        DataStructs.ConvertToNumpyArray(MACCS_val, arr)
        df.loc[len(df)] = arr 
    return df.to_numpy()

# load models and complete prediction for each
outcome_df = pd.DataFrame({'smiles': smiles_list})
current_path = os.path.split(os.path.realpath(__file__))[0]
for index, row in model_meta.iterrows():
    cf = joblib.load(checkpoints_dir + "/classification_model" + row['path'])
    if row['MACCS'] == 1:
        FPs = smiles_to_MACCS(smiles_list)
    else:
        FPs = smiles_to_ECFP(smiles_list, int(row['ECFP_type']), int(row['ECFP_nbits']))
    y_hat_label = cf.predict(FPs)
    y_predict_proba = cf.predict_proba(FPs)
    outcome_df[row['name'] + '_label'] = y_hat_label
    outcome_df[row['name'] + '_prob'] = [i[1] for i in y_predict_proba]

# write output to csv
outcome_df.to_csv(output_file, index=False)
