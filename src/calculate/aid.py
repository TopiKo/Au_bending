'''
Created on 18.11.2014

@author: tohekorh
'''
import os
from ase.io import PickleTrajectory
import numpy as np
#from ase.lattice.surface import fcc111
from ase.md.langevin import Langevin
from ase import units, Atoms
from ase.md import MDLogger
from ase.visualize import view
from ase.lattice.surface import fcc111

def checkAndCreateFolder(path):
    
    if not os.path.exists(path):
        os.makedirs(path)
        
def checkIfExists(T, mdsteps, ind, path,  *args):

    # Check if given trajectory file exists:
    if ind == 'flat':
        path_file     =   path + 'flat/md_data/T=%.0f/md_flat.traj' %T
        #print path_file
    elif ind == 'cylinder':
        R   =   args[0]
        path_file     =   path + 'cylinder/md_data/T=%.0f/md_R%.3f.traj' %(T, R)
    elif ind == 'sphere':
        R   =   args[0]
        path_file     =   path + 'sphere/md_data/T=%.0f/md_R%.3f.traj' %(T, R)
    else: raise

    exists          =   True
    try:
        traj        =   PickleTrajectory(path_file)
        if len(traj) != mdsteps:
            print 'There is differetn number of md steps in file: %s' %path_file 
            exists  =   False
    except IOError:
        print 'There is no file in: %s' %path_file 
        exists      =   False
    return exists

def trim(table):
    
    take        =   np.where(table[:,0] != 0)[0]
    new_table   =   np.zeros((len(take), len(table[0])))
    
    n           =   0
    for i in take:
        new_table[n]    =   table[i]
        n              +=   1
    
    return new_table


def get_Rset(angles, T, params):
    
    # Get the optimized flat strcture, if does not exist -> raise
    atoms       =   get_opmAtomsFlat(T, params)
    L           =   atoms.get_cell().diagonal()
     
    atoms.set_cell(L) 
    
    R           =   np.zeros(len(angles))
    
    # Determine radiuces: R*angle = L[1]
    for i, angle in enumerate(angles):
        R[i]    =   L[1]/angle
    
    return R

def get_AvRadius(traj):
    
    for atoms in traj:
        cell    =   atoms.get_cell()
        R, N    =   0., len(atoms)
        for a in atoms:
            R  +=   np.sqrt(a.position[0]**2 + a.position[1]**2)/N
    
    angle   =   cell[2,0]
    height  =   cell[1,0]
    return R, angle, height

def get_Area(traj):
    
    atoms   =   traj[0]
    cell    =   atoms.get_cell()        
    area    =   np.linalg.norm(np.cross(cell[0], cell[1]))
    return area
    
def get_areaPerAtom(atoms, ind):
    
    if ind == 'flat':
        cell    =  atoms.get_cell()
        area    =  np.linalg.norm(np.cross(cell[0], cell[1]))
        return area/len(atoms)
    
        
    elif ind == 'cylinder':
        cell    =   atoms.get_cell()
        angle   =   cell[2,0]
        height  =   cell[1,0]
        R, N    =   0., len(atoms)
        for a in atoms:
            R  +=   np.sqrt(a.position[0]**2 + a.position[1]**2)/N
        
        area    =   angle*R*height    
        return area/N, R
    else: raise

def get_opmAtomsFlat(T, params):
    
    nm      =   params['mdsteps']
    path    =   params['path']
    if checkIfExists(T, nm, 'flat', path):
        traj   =   PickleTrajectory(path + 'flat/md_data/T=%.0f/md_flat.traj' %T)
        
        return traj[0]

def get_energyPerArea(traj, ind):

    eflatPerArea    =   np.zeros(len(traj))
    R               =   np.zeros(len(traj))
        
    for i, atoms in enumerate(traj):
        e                   =   atoms.get_potential_energy()
        if ind == 'cylinder':
            areaPerAtom, R[i]   =   get_areaPerAtom(atoms, ind)
        else:
            areaPerAtom         =   get_areaPerAtom(atoms, ind)
        eflatPerArea[i]         =   e/(areaPerAtom*len(atoms))
        
        
    return eflatPerArea, np.mean(R) 

def scale_cell_flat(cell_init, s):
    # Scale the unit-cell in x and y - directions
    cell    =   np.zeros((3,3))
    cell[0] =   s*cell_init[0]
    cell[1] =   s*cell_init[1]
    cell[2] =   cell_init[2]
    return cell

def run_dynamics(atoms, opt, params, T, p_trajF, p_logF, ind, *args):
    
    fmax    =   params['fmax']
    dt      =   params['dt']
    fric    =   params['fric']
    mdsteps =   params['mdsteps']
    
    if ind == 'flat':
        
        s                   =   args[0]
        cell_opm            =   args[1]
        init_positions      =   args[2]
        init_cell           =   args[3]
        
        cell                =   scale_cell_flat(cell_opm, s)
        
        # Restore initial state:
        atoms.set_cell(init_cell)
        atoms.positions     =   init_positions
        
        # Set optimal cell
        atoms.set_cell(cell, scale_atoms = True)
        
        # Relax (should be otimal i.e. fmax = 0..)
        opt.run(fmax = fmax, steps = 1000)
        
        # Trajectory for dynamics
        traj        =   PickleTrajectory(p_trajF, 'w', atoms)   
        
        # LAngevin dynamics
        dyn         =   Langevin(atoms, dt*units.fs, units.kB*T, fric)
        
        dyn.attach(MDLogger(dyn, atoms, p_logF, header=True, stress=False,
               peratom=True, mode="w"), interval = 1)
        dyn.attach(traj.write, interval = 1)
        dyn.run(mdsteps)
        traj.close()
        
def get_square_uCell(size):
    
    d           =   2.934
    atoms_flat  =   fcc111('Au',size=(2,2,1),a=np.sqrt(2)*d)
    
     
    idx_sqr     =   [0,2] #[2, 3, 4, 7, 8, 9, 11, 12, 13, 16, 17, 18]
    
    atoms_sqr   =   Atoms()
    for i, a in enumerate(atoms_flat):
        if i in idx_sqr:
            atoms_sqr   +=  a
    
    cell    =   [d, np.sqrt(3)*d, d]   
    
    
    atoms_sqr.center()
    atoms_sqr.translate(-atoms_sqr.positions[0])
    atoms_sqr.set_pbc((True, True, False))
    atoms_sqr.set_cell(cell)
    
    
    
    
    atoms_sqr = atoms_sqr.repeat(size)
    
    
    
    return atoms_sqr
    