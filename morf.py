from os import system, getcwd, makedirs, rename, chdir
from os.path import splitext, isfile
from shutil import copyfile, copytree
import numpy as np
import pandas as pd
import platform
import argparse
import MDAnalysis as mda
from var import *

parser = argparse.ArgumentParser()
parser.add_argument("-gmx", type=str, help="gromacs --- default gmx", default="gmx")
parser.add_argument(
    "-mdrun", type=str, help="mdrun --- default mdrun_tmpi", default="mdrun_tmpi"
)
parser.add_argument("-mother", type=str, help="reference folder to morf")
parser.add_argument("-lig", type=str, help="ligand name .mol2 namefile")
args = parser.parse_args()

system(f"{MATCH} {args.lig}.mol2")
system(
    f"python {charmm2gmx} {getcwd()}/{args.lig}.rtf {getcwd()}/{args.lig}.prm {args.lig}.ff"
)
system(f"obabel -imol2 {args.lig}.mol2 -opdb -O {args.lig}.pdb")
system(
    f"{args.gmx} pdb2gmx --ff {args.lig} -f {args.lig}.pdb -o Ligand_gmx.pdb -p Ligand.top"
)
ligand_ff = f"{args.lig}.ff"
with open(f"{ligand_ff}/ffbonded.itp") as ffbonded, open(
    f"{ligand_ff}/ffnonbonded.itp"
) as ffnonbonded, open(f"Ligand.top") as ligand_top, open(
    f"{args.lig}.itp", "w"
) as ligand_itp:
    ffnonbonded_lines = ffnonbonded.readlines()
    ffbonded_lines = ffbonded.readlines()
    ffnonbonded_lines += ffbonded_lines
    ligand_top_lines = ligand_top.readlines()
    ligand_itp.writelines(ffnonbonded_lines)
    for line in ligand_top_lines:
        if line.startswith(";") or line.startswith("#"):
            pass
        elif "Other" in line:
            line = "UNK   3"
            ligand_itp.write(line)
        elif "[ system ]" in line:
            break
        else:
            ligand_itp.write(line)
    ligand_itp.write(
        """
#ifdef POSRES
#include "../posre.itp"
#endif
"""
    )
print(f"{ligand_itp.name} scritto!")

copytree(args.mother, "madre")  # dovrei creare prima una cartella
chdir("madre")
mother_universe = mda.Universe("step5_input.gro")
ligand_universe = mda.Universe("../Ligand_gmx.pdb")
new_system = mda.Merge(
    mother_universe.select_atoms("all and not resname UNK"), (ligand_universe.atoms)
)
new_system.dimensions = mother_universe.dimensions
new_system.atoms.write("morfed.gro")  # written with dimensions
rename(f"../{args.lig}.itp", "toppar/UNK.itp")
with open("morftopol.top", "w") as topol:
    topol.write(
        """
; Include forcefield parameters
#include "toppar/forcefield.itp"
#include "toppar/UNK.itp"
#include "toppar/PROA.itp"
#include "toppar/PROB.itp"
#include "toppar/POPC.itp"
#include "toppar/POT.itp"
#include "toppar/CLA.itp"
#include "toppar/TIP3.itp"

[ system ]
; Name
Title

[ molecules ]
; Compound	#mols
PROA  	           1
PROB  	           1
POPC  	          40
POT   	           9
CLA   	          13
TIP3  	        3859
UNK   	           1
"""
    )

system("echo q | gmx make_ndx -f morfed.gro -o index.ndx")

new_system.select_atoms("protein or resname UNK").write(
    "index.ndx", name="SOLU", mode="a"
)
new_system.select_atoms("resname TIP3 POT CLA").write(
    "index.ndx", name="SOLV", mode="a"
)
new_system.select_atoms("resname POPC").write("index.ndx", name="MEMB", mode="a")
new_system.select_atoms("protein or (resname UNK) or (resname POPC)").write(
    "index.ndx", name="SOLU_MEMB", mode="a"
)

system("echo 0 | gmx genrestr -f ../Ligand_gmx.pdb -o posre.itp -fc 1000 1000 1000")
system(
    "gmx grompp -f step6.0_minimization.mdp -o step6.0_minimization.tpr -c morfed.gro -r morfed.gro -p morftopol.top -n index.ndx -maxwarn 10"
)
