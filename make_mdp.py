def make_mdp(mdp = 'MD', ns = 5, dt = 0.002, nstxout = 0, nstvout = 0, nstenergy = 0):
    
    if mdp == 'MD':
        with open (f'{mdp}.mdp', 'w') as fout:
            nsteps = int(ns * 1000 / dt)
            save_everyframes = int(nsteps / 1000)
            fout.write(f""";
integrator              = md
dt                      = {dt}
nsteps                  = {nsteps}
; Output control
nstxout     = {nstxout}          
nstvout     = {nstvout}         
nstenergy   = {nstenergy}      
nstlog      = {save_everyframes} 
nstxtcout   = {save_everyframes}  
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
coulombtype             = pme
rcoulomb                = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
;
tcoupl                  = Nose-Hoover
tc_grps                 = SYSTEM
tau_t                   = 1.0
ref_t                   = 303.15
;
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 5.0
compressibility         = 4.5e-5
ref_p                   = 1.0
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = SYSTEM
;
refcoord_scaling        = com
""")

    elif mdp == 'equi':
        with open (f'{mdp}.mdp', 'w') as fout:
            fout.write(f""";
define                  = -DPOSRES
integrator              = md
dt                      = 0.001
nsteps                  = 125000
nstxtcout               = 5000
nstvout                 = 5000
nstfout                 = 5000
nstcalcenergy           = 100
nstenergy               = 1000
nstlog                  = 1000
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
coulombtype             = pme
rcoulomb                = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
;
tcoupl                  = Nose-Hoover
tc_grps                 = SYSTEM
tau_t                   = 1.0
ref_t                   = 303.15
;
constraints             = h-bonds
constraint_algorithm    = LINCS
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = SYSTEM
;
gen-vel                 = yes
gen-temp                = 303.15
gen-seed                = -1
;
refcoord_scaling        = com""")
 
    elif mdp == 'mini':
        with open (f'{mdp}.mdp', 'w') as fout:
            fout.write(f""";
define                  = -DPOSRES
integrator              = steep
emtol                   = 1000.0
nsteps                  = 5000
nstlist                 = 10
cutoff-scheme           = Verlet
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = pme
rcoulomb                = 1.2
;
constraints             = h-bonds
constraint_algorithm    = LINCS""")

    elif mdp == 'ions':
        with open (f'{mdp}.mdp', 'w') as fout:
            fout.write(f"""
; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 200           ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions""")
    else: print('something wrong?')
    return fout.name
