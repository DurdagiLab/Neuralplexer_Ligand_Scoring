#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 11:22:59 2024

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

import os
import glob
import shutil
import multiprocessing

base_folder_path = '/path for NeuralPlexer output/'
os.chdir(base_folder_path)

def process_folder(folder_path):
    os.chdir(folder_path)
    files= glob.glob('*')
    for folder_index in files:
        f_path = os.path.join(folder_path, f'{folder_index}')
        if os.path.isdir(f_path):
            print(f"Processing folder: {f_path}")
            os.chdir(f_path)
            os.system("$SCHRODINGER18/utilities/sdconvert -isd lig_rank1.sdf -omae lig_rank1.mae")
            os.system("$SCHRODINGER18/utilities/pdbconvert -ipdb protein_rank1.pdb -omae protein_rank1.mae")
            os.system("cat protein_rank1.mae lig_rank1.mae > complex_pv.mae")
            os.system("$SCHRODINGER18/run pv_convert.py -mode merge complex_pv.mae")
            os.system(f"$SCHRODINGER18/utilities/prepwizard -disulfides -fillsidechains -propka_pH 7.4 -epik_pH 7.4 -epik_pHt 0.5 -noimpref -j {folder_index} complex-out_complex.mae prep-complex.mae -WAIT")
            os.system(f'$SCHRODINGER18/utilities/generate_glide_grids -rec_file prep-complex.mae -lig_asl "res.ptype UNK" -j {folder_index} -SAVE -WAIT')
            os.system("$SCHRODINGER18/run pv_convert.py prep-complex.mae -mode split_ligand -s")
            with open('glide.in', 'w') as glide:
                glide.write(f'CALC_INPUT_RMS   True\nDOCKING_METHOD   mininplace\nFORCEPLANAR   True\nGRIDFILE   {f_path}/{folder_index}-gridgen.zip\nLIG_MAECHARGES   True\nLIGANDFILE   {f_path}/prep-complex-out1_lig.mae\nPOSTDOCK   False\nPRECISION   SP\nREWARD_INTRA_HBONDS   True')
            os.system('$SCHRODINGER18/glide glide.in -OVERWRITE -adjust -HOST "localhost:24" -TMPLAUNCHDIR -ATTACHED -WAIT')
            os.system(f'$SCHRODINGER18/utilities/proplister -p s_i_glide_gridfile r_i_docking_score r_i_glide_rmsd_to_input -c -u glide_pv.maegz -o {folder_index}-score.txt')
        else:
            print(f"Folder not found: {f_path}")


def main():
    process_folder(base_folder_path+'subfolder_1/')

if __name__ == "__main__":
    main()
