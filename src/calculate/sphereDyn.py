
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





def trajToExtTraj(traj, n, *args):
    
    R, T, angle1, angle2, n1, n2, path = args
   
    path_traj   =   path + 'sphere/md_data/T=%.0f/' %T
    checkAndCreateFolder(path_traj)
    traj_ext    =   PickleTrajectory(path_traj + 'md_Ext_R%.3f.traj' %R, 'w') 
    
    for a in traj:
        c_symb  =   a.get_chemical_symbols()
        posits  =   a.positions
        atoms   =   Atoms(c_symb, positions = posits, container='Sphere')
        atoms.set_container(angle1=angle1,angle2=angle2,n1=n1,n2=n2,mode=4)
        traj_ext.write(atoms.extended_copy(n))


def sphereDynamics(R, T, params):
    
    d       =   params['d']         #2.934
    nk      =   params['nk']        #20
    fmax    =   params['fmax']      #1E-2
    fric    =   params['fric']      #0.002        # equilibration time tau = 10 fs/fric (= mdsteps*dt) 
    dt      =   params['dt']        #5.0
    mdsteps =   params['mdsteps']   #int(10/fric/dt)
    path    =   params['path'] #'/space/tohekorh/Au_bend/files/'
    
    print 'radius',R
    name    =   '%.1f' %R
    
    # ATOMS
    atoms   =   fcc111('Au',size=(2,2,1),a=np.sqrt(2)*d)
    length1 =   np.linalg.norm( atoms.get_cell()[0] )
    length2 =   np.linalg.norm( atoms.get_cell()[1] )
    n1      =   (0,1,0)
    n2      =   (-np.cos(np.pi/6),np.sin(np.pi/6),0)
    angle1  =   length1/R
    angle2  =   length2/R
    
    atoms.translate((0,0,R))
    
    atoms   =   Atoms(atoms=atoms,container='Sphere')
    atoms.set_container(angle1=angle1,angle2=angle2,n1=n1,n2=n2,mode=4)
    
    # FOLDERS
    path_md =   path + 'sphere/md_data/T=%.0f/' %T
    path_opm=   path + 'sphere/opm/T=%.0f/' %T

    checkAndCreateFolder(path_md)
    checkAndCreateFolder(path_opm)
    
    # CALCULATOR
    calc    =   Hotbit(SCC=False, kpts=(nk,nk,1),txt= path_opm + 'optimization_%s.cal' %name,)
    atoms.set_calculator(calc)
    
    # RELAX
    opt     =   BFGS(atoms, trajectory= path_opm + 'optimization_%s.traj' %name)
    opt.run(fmax=fmax,steps=300)
    
    
    # DYNAMICS
    traj        =   PickleTrajectory(path_md + 'md_R%.3f.traj' %R, 'w', atoms)   
    dyn         =   Langevin(atoms, dt*units.fs, units.kB*T,fric)
    
    dyn.attach(MDLogger(dyn, atoms, path_md + 'md_R%.3f.log' %R, header=True, stress=False,
           peratom=True, mode="w"), interval = 1)
    dyn.attach(traj.write)
    dyn.run(mdsteps)
    traj.close()
    
    # load the dynamics back..
    traj        =   PickleTrajectory(path_md + 'md_R%.3f.traj' %R)
    trajToExtTraj(traj, (2, 2, 1), R, T, angle1, angle2, n1, n2, path)
    
    #cohesion    =   np.mean([atoms.get_potential_energy() / len(atoms) for atoms in traj]) 
    
    #cohesion    =   atoms.get_potential_energy()/len(atoms)
    #atoms2      =   atoms.extended_copy((2,2,1))
    #view(atoms2)
    
    #return cohesion
    
    
    
    
    