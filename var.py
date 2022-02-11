#var.py

MATCH = '$HOME/htmd/MATCH_RELEASE/MATCH/scripts/MATCH.pl'
charmm2gmx = '$HOME/htmd/charmm2gromacs-pvm3.py'
null = '2> /dev/null'
columns_name = ['ligand','status', 'contacts_average', 'all_contacts_plot','RMSD', 'coordination_plot', 'smiles', 'made by']
colors = {'Hydrophobic': 'purple', 'Pi_stacking_T': 'black', 'Pi_stacking_P': 'coral', 'Salt_Bridge': 'mediumblue', 'H_bond':'aqua', 'Pi_Cation':'firebrick', 'Water_Bridge':'gold'}

METAL_IONS = ['CA', 'CO', 'MG', 'MN', 'FE', 'CU', 'ZN', 'FE2', 'FE3', 'FE4', 'LI', 'RB', 'SR', 'CS', 'BA',
              'CR', 'NI', 'FE1', 'NI', 'RU', 'RU1', 'RH', 'RH1', 'PD', 'AG', 'CD', 'LA', 'W', 'W1', 'OS', 'IR', 'PT',
              'PT1', 'AU', 'HG', 'CE', 'PR', 'SM', 'EU', 'GD', 'TB', 'YB', 'LU', 'AL', 'GA', 'IN', 'SB', 'TL', 'PB', 'MN2P']
