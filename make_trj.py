
import MDAnalysis as mda
from MDAnalysis.analysis import align
import argparse

def _parse_CLAs() -> argparse.ArgumentParser:
    """
    Parser Command Line Arguments

    Returns:
        argparse.ArgumentParser : parsed command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-topol", 
        type=str, 
        help="file pdb,psf,tpr --- default MD.pdb",
        default="MD.pdb"
    )
    parser.add_argument(
        "-trj", 
        type=str, 
        help="file xtc,trr,dcd --- default MD.xtc", 
        default="MD.xtc"
    )
    parser.add_argument(
        "-step", 
        type=int, 
        help="value of step --- default 10", 
        default="10"
    )
    parser.add_argument(
        "-alignto", 
        type=str, 
        help="select atoms to align", 
        default="name CA"
    )
    parser.add_argument(
        "-filename", 
        type=str, 
        help="output trajectory", 
        default="aligned.dcd"
    )
    parser.add_argument(
        "-outgroup", 
        type=str, 
        help="select outgroup", 
        default="protein"
    )
    return parser

def myalign(args: argparse.ArgumentParser):
    """
    Align trajectory to reference, and export output selection"""
    reference = mda.Universe(args.topol)
    specimen = mda.Universe(args.topol, args.trj)

    align.AlignTraj(specimen,
                    reference,
                    select = args.alignto,
                    filename = "temp." + args.filename[-3:],
                    match_atoms = True,
                    ).run(step = args.step)
    
    u = mda.Universe(args.topol, args.trj)
    outgroup = u.select_atoms(args.outgroup)
    with mda.Writer("temp." + args.filename[-3:], outgroup.n_atoms) as W:
        for ts in u.trajectory:
            W.write(args.filename)
    return args.filename


if __name__ == "__main__":
    parser = _parse_CLAs()
    args = parser.parse_args()

    myalign(args)