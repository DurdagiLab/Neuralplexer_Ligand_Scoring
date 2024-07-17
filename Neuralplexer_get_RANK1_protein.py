#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 13:26:55 2024

@author: Ehsan.Sayyah
"""
import numpy as np
from rdkit import Chem
import os
import argparse
import glob
import shutil
import time
import subprocess

def read_sdf(file_path):
    # Read SDF file and return a list of molecules
    return [mol for mol in Chem.SDMolSupplier(file_path) if mol is not None]

def extract_model(pdb_file_path, match_index):
    coordinates = []
    with open(pdb_file_path, 'r') as file:
        model_found = False
        for line in file:
            if line.startswith('MODEL'):
                current_model = int(line.split()[1])
                if current_model == match_index:
                    model_found = True
                else:
                    model_found = False
            if model_found==True:
                if line.startswith('ATOM') or line.startswith('ENDMDL') or line.startswith('TER'):
                    coordinates.append(line)
    return coordinates

def get_coordinates(molecule):
    coords = []
    conf = molecule.GetConformer()
    for atom in molecule.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        coords.append([pos.x, pos.y, pos.z])
    return np.array(coords)

def compare_coordinates(coords1, coords2, tolerance=1e-3):
    return np.allclose(coords1, coords2, atol=tolerance)

def process_folder(folder_path, folder_index):
    # Paths to files in the current folder
    lig_all_path = os.path.join(folder_path, 'lig_all.sdf')
    lig_rank1_path = os.path.join(folder_path, 'lig_rank1.sdf')
    pdb_file_path = os.path.join(folder_path, 'prot_all.pdb')
    
    # Read ligands
    lig_all = read_sdf(lig_all_path)
    lig_rank1 = read_sdf(lig_rank1_path)
    
    for i, lig1 in enumerate(lig_rank1):
        # Get coordinates of the ligand in lig_rank1
        coords1 = get_coordinates(lig1)
        
        # Find the matching ligand index in lig_all
        match_index = -1
        for j, lig2 in enumerate(lig_all):
            coords2 = get_coordinates(lig2)
            if compare_coordinates(coords1, coords2):
                match_index = j + 1  # Model indices are 1-based
                break
        
        if match_index == -1:
            print(f"No match found for ligand {i} in {lig_rank1_path}")
            continue
        
        # Extract the CA coordinates for the matched model
        ca_coords = extract_model(pdb_file_path, match_index)
        pdb_coord= ''.join(ca_coords)
        with open(folder_index+'/protein_rank1.pdb', 'w') as pdb:
            pdb.write(pdb_coord)

def main():
    # Base path to folders
    base_folder_path = '/media/arma/DATA/S-Ehsan/NeuralPLexer/HIV_plexer/'
    os.chdir(base_folder_path)
    files= glob.glob('line*')
    # Initialize the final matrix
    # num_folders = 140
    for folder_index in files:
        folder_path = os.path.join(base_folder_path, f'{folder_index}')
        if os.path.isdir(folder_path):
            print(f"Processing folder: {folder_path}")
            process_folder(folder_path, folder_index)
        else:
            print(f"Folder not found: {folder_path}")
    

if __name__ == "__main__":
    main()
