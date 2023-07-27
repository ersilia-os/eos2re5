#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tues Jul 26 2022

@author: Samuel Volk
"""
# point to admetlab/chemopy paths
import os, sys
path_root = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..")
sys.path.append(os.path.abspath(path_root))

# additional imports
from sklearn.externals import joblib
import ast
import numpy as np
import pandas as pd
from chemopy.src.pychem import constitution, topology, moran, moe, \
    molproperty, connectivity, bcut, estate, basak, moreaubroto
from chemopy.src.pychem.pychem import PyChem2d
from rdkit.Chem import DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit import Chem

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
file_path = os.path.dirname(os.path.abspath(__file__))

# hide gpus
os.environ["CUDA_VISIBLE_DEVICES"]="-1"

# configure relevant resource paths
checkpoints_dir = os.path.abspath(os.path.join(file_path, "..", "..", "checkpoints"))
framework_dir = os.path.abspath(os.path.join(file_path, ".."))

# read smiles input and model metadata csv files
smiles_df = pd.read_csv(input_file)
smiles_list = smiles_df[smiles_df.columns[0]].tolist()
model_meta = pd.read_csv(os.path.join(framework_dir, 'model_metadata.csv'))

# function to parse descriptor list needed for some models and assemble input
def smiles_to_descriptors(smiles_list, descriptor_list):
    df = pd.DataFrame(columns=descriptor_list)
    pc2D = PyChem2d()
    for i in smiles_list:
        mol_rdkit = Chem.MolFromSmiles(i)
        pc2D.ReadMolFromSmile(i)
        des = constitution.GetConstitutional(mol_rdkit)
        des.update(moe.GetMOE(mol_rdkit)) # concat new descriptors
        des.update(topology.GetTopology(mol_rdkit))
        des.update(connectivity.GetConnectivity(mol_rdkit))
        des.update(moran.GetMoranAuto(mol_rdkit))
        des.update(moreaubroto.GetMoreauBrotoAuto(mol_rdkit))
        des.update(molproperty.GetMolecularProperty(mol_rdkit))
        des.update(estate.GetEstate(mol_rdkit))
        des.update(bcut.GetBurden(mol_rdkit))
        des.update(basak.Getbasak(mol_rdkit))
        des.update(pc2D.GetCharge())
        des.update(pc2D.GetKappa())
        df.loc[len(df)] = [des[j] for j in descriptor_list] 
    return df.to_numpy()   
 
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
outcome_df = pd.DataFrame()
is_not_ECFP = model_meta['ECFP_type'].isnull().to_numpy()
is_not_MACCS = model_meta['MACCS'].isnull().to_numpy()
for index, row in model_meta.iterrows():
    # load model checkpoint
    cf = joblib.load(checkpoints_dir + row['path'])

    # form input based on model input type
    if is_not_MACCS[index] and is_not_ECFP[index]: # model takes chemical descriptors as input
        descriptor_list = ast.literal_eval(row['descriptors'])
        input = smiles_to_descriptors(smiles_list, descriptor_list)
       
    elif is_not_MACCS[index]: # model takes ECFP as input
        input = smiles_to_ECFP(smiles_list, int(row['ECFP_type']), int(row['ECFP_nbits']))
    
    else: # model takes MACCS as input
        input = smiles_to_MACCS(smiles_list)

    # generate predictions based on model type (regression/classification)
    if row['model_type'] == "classification": # classification model
        y_hat = cf.predict(input)
        y_proba = cf.predict_proba(input)
        outcome_df[row['name'] + '_label'] = y_hat
        outcome_df[row['name'] + '_prob'] = [i[1] for i in y_proba]
    else: # regression model
        y_pred = cf.predict(input)
        outcome_df[row['name'] + '_pred'] = y_pred

# write output to csv
outcome_df.to_csv(output_file, index=False)
