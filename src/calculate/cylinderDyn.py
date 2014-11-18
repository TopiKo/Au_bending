
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
from ase.visualize import view
#mdsteps =   10





def trajToExtTraj(traj, n, *args):
    
    R, T, angle, H, path = args
   
    path_traj   =   path + 'cylinder/md_data/T=%.0f/' %T
    checkAndCreateFolder(path_traj)
    
    traj_ext = PickleTrajectory(path_traj + 'md_Ext_R%.3f.traj' %R,'w') 
    for a in traj:
        c_symb  =   a.get_chemical_symbols()
        posits  =   a.positions
        atoms   =   Atoms(c_symb, positions = posits, container='Wedge')
        atoms.set_container(angle=angle, height=H, physical=False, pbcz=True)
        traj_ext.write(atoms.extended_copy(n))


def cylinderDynamics(R, T, params):
    
    d       =   params['d']         #2.934
    nk      =   params['nk']        #20
    fmax    =   params['fmax']      #1E-2
    fric    =   params['fric']      #0.002        # equilibration time tau = 10 fs/fric (= mdsteps*dt) 
    dt      =   params['dt']        #5.0
    mdsteps =   params['mdsteps']   #int(10/fric/dt)
    path    =   params['path'] #'/space/tohekorh/Au_bend/files/'
    
    print 'radius',R
    name    =   '%.1f' %R
    
    #
    
    atoms       =   fcc111('Au',size=(2,2,1),a=np.sqrt(2)*d)
    L           =   atoms.get_cell().diagonal()
    
    atoms.set_cell(L) 
    
    angle       =   L[1]/R
    atoms.rotate('y', np.pi/2)
    atoms.translate((-atoms[0].x, 0, 0) )
    atoms.translate((R,0,0))
    atoms       =   Atoms(atoms=atoms,container='Wedge')
    atoms.set_container(angle=angle,height=L[0],physical=False,pbcz=True)
    #view(atoms.extended_copy((8,1,2)))
    
    # FOLDERS
    path_opm    =   path + 'cylinder/opm/T=%.0f/' %T
    path_md     =   path + 'cylinder/md_data/T=%.0f/' %T
    
    checkAndCreateFolder(path_opm)
    checkAndCreateFolder(path_md)
    
    # proper number of kappa points in angle
    m           =   int( round(2*np.pi/angle) )
    # CALCULATOR
    calc        =   Hotbit(SCC=False, kpts=(m,1,nk), \
                    txt= path_opm + 'optimization_%s.cal' %name)
    atoms.set_calculator(calc)
    
    opt         =   BFGS(atoms,trajectory= path_opm + 'optimization_%s.traj' %name)
    opt.run(fmax=fmax,steps=300)
    
    
    # DYNAMICS
    traj        =   PickleTrajectory(path_md + 'md_R%.3f.traj' %R,'w',atoms)   
    
    dyn         =   Langevin(atoms, dt*units.fs, units.kB*T,fric)
    dyn.attach(MDLogger(dyn, atoms, path_md + 'md_R%.3f.log' %R, header=True, stress=False,
           peratom=True, mode="w"), interval = 1)
    
    dyn.attach(traj.write)
    dyn.run(mdsteps)
    traj.close()
    
    # load the dynamics back..
    traj        =   PickleTrajectory(path_md + 'md_R%.3f.traj' %R)
    trajToExtTraj(traj, (4, 1, 2), R, T, angle, L[0], path)
    
    #cohesion    =   np.mean([atoms.get_potential_energy() / len(atoms) for atoms in traj]) 
    
    #cohesion    =   atoms.get_potential_energy()/len(atoms)
    #atoms2      =   atoms.extended_copy((2,2,1))
    #view(atoms2)
    
    #return cohesion
    
    
    
    
    