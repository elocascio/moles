import os
from ex.DNAze_pdb import DNAze_pdb

os.system(f"gmx editconf -f MD.tpr -o reference.pdb")
DNAze_pdb("./reference.pdb")
print("finito")
