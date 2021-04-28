#var.py

gmx = 'gmx'
MATCH = '$HOME/htmd/MATCH_RELEASE/MATCH/scripts/MATCH.pl'
charmm2gmx = '$HOME/htmd/charmm2gromacs-pvm3.py'
mdrun = 'mdrun_tmpi'
name = 'MG-REX'
null = '2> /dev/null'
columns_name = ['ligand', 'contacts_average', 'all_contacts', 'coordination_plot','made by']
colors = {'Hydrophobic': 'purple', 'Pi_stacking_T': 'black', 'Pi_stacking_P': 'coral', 'Salt_Bridge': 'mediumblue', 'H_bond':'aqua', 'Pi_Cation':'firebrick', 'Water_Bridge':'gold'}
