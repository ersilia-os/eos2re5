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
import os

smiles_list = [
    "CC1(C)OC2CC(=O)OCC23C1CC(=O)C1(C)C3CCC2(C)C(c3ccoc3)OC(=O)C3OC321",
    "OCC1OC(OC2(CBr)OC(CBr)C(O)C2O)C(O)C(O)C1Br",
    "CC1C2CC(C)CCN2C2CC3C4CC=C5CC(OC6OC(CO)C(O)C(OC7OC(CO)C(O)C(O)C7O)C6OC6OC(C)C(O)C(O)C6O)CCC5(C)C4CCC3(C)C12",
    "COc1ccc2c(c1)S(=O)(=O)NC2=O",
    "CCOc1ccc(CCC(=O)c2c(O)cccc2O)cc1O",
    "COc1ccc(NCc2ccccc2)cc1O",
    "COc1cc2c(cc1OC)C13CCN4CC5=CCOC6CC(=O)N2C1C6C5CC43",
    "NC(Cc1c[nH]c2cc(Cl)ccc12)C(=O)O",
    "COc1ccc(C2Cc3cccc(OC(C)=O)c3C(=O)O2)cc1OC(C)=O"
]

model_meta = pd.read_csv('model_metadata.csv')

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
    cf = joblib.load(current_path+'/classification_model'+row['path'])
    if row['MACCS'] == 1:
        FPs = smiles_to_MACCS(smiles_list)
    else:
        FPs = smiles_to_ECFP(smiles_list, int(row['ECFP_type']), int(row['ECFP_nbits']))
    y_hat_label = cf.predict(FPs)
    y_predict_proba = cf.predict_proba(FPs)
    outcome_df[row['name'] + '_label'] = y_hat_label
    outcome_df[row['name'] + '_prob'] = [i[1] for i in y_predict_proba]

# write output to csv
outcome_df.to_csv("outcome.csv", index=False)