import pandas as pd
from os import chdir, makedirs
from shutil import copytree
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "-table", type=str, help="data table, tab spaced")
parser.add_argument("-b", "-best", type=int, help="best ligand to take", default = 100)
args = parser.parse_args()

makedirs(f'BEST_{args.best}', exist_ok=True)

df = pd.read_table({args.table})
df = df.head(args.best)

for element in df['ligand'].values:
    copytree(element, f'BEST_{args.best}/{element}')
