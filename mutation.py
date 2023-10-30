from pymol import cmd
from itertools import cycle
import argparse
from Misc.moles import init


def mutation(pdb, residues, mutation, pdb_out):
    cmd.wizard("mutagenesis")
    cmd.load(pdb)
    cmd.refresh_wizard()

    res_list = residues.split(",")
    res_cycle = cycle(res_list)
    aminoacid_list = mutation.split(",")
    amino_cycle = cycle(aminoacid_list)

    if len(res_list) != len(aminoacid_list):
        print("unmatchable items!")
        exit()

    for n in range(len(res_list)):
        chain, res = next(res_cycle).split("-")
        cmd.get_wizard().do_select(f"{chain}/{res}/")
        cmd.get_wizard().set_mode(next(amino_cycle).upper())
        cmd.get_wizard().apply()

    cmd.set_wizard()
    cmd.save(pdb_out)
    return pdb_out


if __name__ == "__main__":

    init()
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, help="input pdb file")
    parser.add_argument(
        "-r", type=str, help='residue to mutate - syntax "A-588,B-577,C-23"'
    )
    parser.add_argument(
        "-a",
        type=str,
        help='aminoacid to substitute, 3 letter, comma separate. IN ORDER - ex. "ANS,PHE,ARG,TRP"',
    )
    parser.add_argument("-o", type=str, help="output pdb file name")
    args = parser.parse_args()

    mutation(args.i, args.r, args.a, args.o)
