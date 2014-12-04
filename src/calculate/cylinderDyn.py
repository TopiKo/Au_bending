
from ase import *
from hotbit import *
import numpy as np
from ase.optimize import BFGS
from ase.io import PickleTrajectory
from ase.md.langevin import Langevin
from ase.md import MDLogger
from aid import checkAndCreateFolder, get_opmAtomsFlat
from ase.visualize import view






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
    
    nk      =   params['nk']        #20
    fmax    =   params['fmax']      #1E-2
    fric    =   params['fric']      #0.002         
    dt      =   params['dt']        #5.0
    mdsteps =   params['mdsteps']   #int(10/fric/dt)
    path    =   params['path']      #'/space/tohekorh/Au_bend/files/'
    
    print 'radius',R, '\n'
    name    =   '%.1f' %R
    
    # Get the optimal flat atoms , unit-cell config for this temperature:
    atoms       =   get_opmAtomsFlat(T, params)

    # Take the diagonal of the cell:
    L           =   atoms.get_cell().diagonal()

     
    atoms.set_cell(L) 
    
    # The wedge angle is:
    angle       =   L[1]/R
    
    # fiddle with atoms:
    atoms.rotate('y', np.pi/2)
    atoms.translate((-atoms[0].x, 0, 0) )
    
    #view(atoms)
    
    # Initial map for atoms, to get the on surface of cylinder:
    phi_max     =   angle
    for a in atoms:
        r0      =   a.position
        phi     =   r0[1]/L[1]*phi_max
        a.position[0]   =   R*np.cos(phi)
        a.position[1]   =   R*np.sin(phi)

    #view(atoms)
    
    # proper number of kappa points in angle direction, if angle >= 2*pi/3 the use genuine 
    # physical boundary conditions and proper kappa-points. With 64 atoms in unit cell souhld not 
    # happen as the radius -> very small if angle = 2*pi/3 or larger.
    
    # Check that the unit-cell angle is 2*pi/integer. This can be removed! 
    if (2*np.pi/angle)%1 > 0:
        raise 
    
    # This does not have any effect unless angle >= 2*pi/3:
    m           =   int(round(2*np.pi/angle))
    physical    =   False
    if m <= 3:
        nk1 =   m
        physical = True
    else:
        nk1     =   nk
    
    # Set up the wedge container:
    atoms       =   Atoms(atoms = atoms, container = 'Wedge')
    atoms.set_container(angle = angle, height = L[0], physical = physical, pbcz = True)
    
    # Check that everything looks good:
    #view(atoms.extended_copy((2,1,2)))
    
    
    # FOLDERS
    path_opm    =   path + 'cylinder/opm/T=%.0f/' %T
    path_md     =   path + 'cylinder/md_data/T=%.0f/' %T
    
    checkAndCreateFolder(path_opm)
    checkAndCreateFolder(path_md)
    
    
    # CALCULATOR
    calc        =   Hotbit(SCC=False, kpts=(nk1, 1, nk), physical_k = physical, \
                    txt= path_opm + 'optimization_%s.cal' %name)
    atoms.set_calculator(calc)
    
    # RELAX
    opt         =   BFGS(atoms, trajectory= path_opm + 'optimization_%s.traj' %name)
    opt.run(fmax=fmax, steps=1000)
    
    
    # DYNAMICS
    traj        =   PickleTrajectory(path_md + 'md_R%.3f.traj' %R,'w',atoms)   
    
    # Perform dynamics
    dyn         =   Langevin(atoms, dt*units.fs, units.kB*T, fric)
    dyn.attach(MDLogger(dyn, atoms, path_md + 'md_R%.3f.log' %R, header=True, stress=False,
           peratom=True, mode="w"), interval = 1)
    
    dyn.attach(traj.write)
    dyn.run(mdsteps)
    traj.close()
    
    # load the dynamics back, to make exteded trajectory:
    traj        =   PickleTrajectory(path_md + 'md_R%.3f.traj' %R)
    trajToExtTraj(traj, (2, 1, 2), R, T, angle, L[0], path)
    
    
    
    
    
    