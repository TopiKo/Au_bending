
from ase import *
from hotbit import *
import numpy as np
from ase.lattice.surface import fcc111
from ase.optimize import BFGS
#from ase import io
#from ase.visualize import view
from ase.io import PickleTrajectory
#from ase import units
from ase.md.langevin import Langevin
from ase.md import MDLogger
from aid import checkAndCreateFolder


#mdsteps =   10


def flatDynamics(T, params):
    
    d       =   params['d']         #2.934
    nk      =   params['nk']        #20
    fmax    =   params['fmax']      #1E-2
    fric    =   params['fric']      #0.002        # equilibration time tau = 10 fs/fric (= mdsteps*dt) 
    dt      =   params['dt']        #5.0
    mdsteps =   params['mdsteps']   #int(10/fric/dt)
    path    =   params['path']      #'/space/tohekorh/Au_bend/files/'
    
    
    atoms       =   fcc111('Au',size=(2,2,1),a=np.sqrt(2)*d)
    
    # FOLDERS
    path_opm    =   path + 'flat/opm/T=%.0f/' %T
    path_md     =   path + 'flat/md_data/T=%.0f/' %T
    
    checkAndCreateFolder(path_md)
    checkAndCreateFolder(path_opm)
    
    # CALCULATOR
    calc        =   Hotbit(SCC=False, kpts=(nk,nk,1),  \
                           txt= path_opm + 'optimization_flat.cal')
    atoms.set_calculator(calc)
    
    # RELAX
    opt         =   BFGS(atoms, trajectory= path_opm + 'optimization_flat.traj')
    opt.run(fmax=fmax,steps=300)
    
    
    # DYNAMICS
    traj        =   PickleTrajectory(path_md + 'md_flat.traj', 'w', atoms)   
    dyn         =   Langevin(atoms, dt*units.fs, units.kB*T,fric)
    
    dyn.attach(MDLogger(dyn, atoms, path_md + 'md_flat.log', header=True, stress=False,
           peratom=True, mode="w"), interval = 1)
    dyn.attach(traj.write)
    dyn.run(mdsteps)
    traj.close()
    
    # load the dynamics back..
    #traj        =   PickleTrajectory(path_md)
    
    #cohesion    =   np.mean([atoms.get_potential_energy() / len(atoms) for atoms in traj]) 
    #trajToExtTraj(traj, (2, 2, 1), T, path)
    
    #cohesion    =   atoms.get_potential_energy()/len(atoms)
    #atoms2      =   atoms.extended_copy((2,2,1))
    #view(atoms2)
    
    #return cohesion
    
    
    
    
    