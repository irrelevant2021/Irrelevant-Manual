#!/bin/env python3
#
# This is an example how to use the ExternalProcess scoring component using
# QSARtuna, see DOI 10.1021/acs.jcim.4c00457.  The scripts expects a list of
# SMILES from stdin and will # write a JSON string to stdout.
#
# QSARtuna code at https://github.com/MolecularAI/QSARtuna.
#
# [[component.ExternalProcess.endpoint]]
# name = "QSARtuna model"
# weight = 0.6
#
# # Run Qptuna in its own environment
# # The --no-capture-output is necessary to pass through stdout from REINVENT4
# # In multi endpoint scenarios replace executable with "/dev/null" in all
# # endpoints except the first
# params.executable = "/home/user/miniconda3/condabin/mamba"
# params.args = "run --no-capture-output -n qsartuna /path/to/run-qsartuna.py model_filename
# params.property = "clint"
#
# # Don't forget the transform if needed!
#

#shapescreen toml
#[[stage.scoring.component]]
#[stage.scoring.component.ExternalProcess]
#[[stage.scoring.component.ExternalProcess.endpoint]]
#name = "Shape Screen"
#weight = 1.0
#params.executable = "/home/xfusion/miniforge3/condabin/mamba"
##sys.argv[0]是run-shapescreen.py; sys.argv[1]是ref.sdf文件,不带.sdf
#params.args = "run --no-capture-output -p /home/xfusion/miniforge3/envs/reinvent4 /home/xfusion/REINVENT4/support/run-shapescreen.py DY-0028"

import sys
import json
import pickle

import pandas as pd
from rdkit import Chem
import os
import subprocess

import numpy as np

smiles = [smiles.strip() for smiles in sys.stdin]
#smiles = ['c1ccccc1'] #测试

# Everything from here to END is specific to the scorer

ref_sdf = sys.argv[1] + '.sdf'
shape_maegz = sys.argv[1] + '_align.maegz'
shape_sdf = sys.argv[1] + '_align.sdf'
shape_csv = sys.argv[1] + '_align.csv'

writer = Chem.SDWriter('temp.sdf') #将来自stdin的smiles转为sdf,每个step的临时temp.sdf会被覆盖
for smile in smiles:
    mol = Chem.MolFromSmiles(smile)
    if mol is None:
        mol = Chem.MolFromSmiles("CCC") #无效的smiles，被替换为占位符,丙烷,保持stdin和stdout数量一致
    writer.write(mol)
writer.close()

schrodinger_path = os.environ.get('SCHRODINGER') #获取schrodinger路径

executable_ligprep = os.path.join(schrodinger_path, 'ligprep')
command_ligprep = [
        executable_ligprep,
        '-isd', 'temp.sdf',
        '-osd', 'ligprep.sdf',
        '-g',
        '-s', '1',
        '-i', '1',
        '-HOST', 'localhost:128',
        '-NJOBS', '50',
        '-WAIT'
    ]
#DEVNULL非常重要,reinvnet通过获取输出来传递数据,所以除了json不能有其他输出
subprocess.run(command_ligprep, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

executable_shape = os.path.join(schrodinger_path, 'shape_screen')
command_shape = [
        executable_shape,
        '-shape', ref_sdf,
        '-screen', 'ligprep.sdf',
        '-flex',
        '-max', '50',
        '-norm', '1',
        '-HOST', 'localhost:128',
        '-TMPLAUNCHDIR',
        '-WAIT'
    ]
subprocess.run(command_shape, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

executable_convert = os.path.join(schrodinger_path, 'utilities/structconvert')
command_sdf = [
        executable_convert,
        shape_maegz,
        shape_sdf
        ]
command_csv = [
        executable_convert,
        shape_sdf,
        shape_csv
        ]
subprocess.run(command_sdf, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
subprocess.run(command_csv, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

df = pd.read_csv(shape_csv)
df_scores = df['r_phase_Shape_Sim']
scores = []
for score in df_scores:
    scores.append(float(score))

# END


# Format the JSON string for REINVENT4 and write it to stdout
# Replace "predictions" if wanted and set params.property accordingly
# Multiple endpoints possible
data = {"version": 1, "payload": {"predictions": list(scores)}}

print(json.dumps(data))
