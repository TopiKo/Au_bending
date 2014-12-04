
from ase import *
from hotbit import *
import numpy as np
from ase.optimize import BFGS
from ase.visualize import view
from ase.io import PickleTrajectory
from aid import checkAndCreateFolder, get_square_uCell, scale_cell_flat, run_dynamics 
from scipy.optimize import fmin, curve_fit
import matplotlib.pyplot as plt





def flatDynamics(T, params):
    
    # Parameters are read from the parameters...
    nk          =   params['nk']        #20
    fmax        =   params['fmax']      #1E-2
    path        =   params['path']      #'/space/tohekorh/Au_bend/files/'
    uzsize      =   params['uz_size']
    therm_step  =   params['thermN']    # from this onward system is thermalized
    
    # Get square fcc unit cell 64 atoms = (8,4,1)
    atoms           =   get_square_uCell(uzsize)
    
    # Store the initial state
    init_positions  =   atoms.positions.copy()
    init_cell       =   atoms.get_cell().copy()
    
    # Create necessary folders 
    path_opm    =   path + 'flat/opm/T=%.0f/' %T
    path_md     =   path + 'flat/md_data/T=%.0f/' %T
    
    checkAndCreateFolder(path_md)
    checkAndCreateFolder(path_opm)
    
    # CALCULATOR
    calc        =   Hotbit(SCC=False, kpts=(nk,nk,1),  \
                           txt= path_opm + 'optimization_flat.cal')
    atoms.set_calculator(calc)
    
    # BFGS object
    opt         =   BFGS(atoms, trajectory= path_opm + 'optimization_flat.traj')
    
    
    # Initial optimization of the unit cell size
    def optim_cell(s):
        
        cell        =   scale_cell_flat(init_cell, s)
        atoms.set_cell(cell, scale_atoms = True)
        opt.run(fmax)
        
        return atoms.get_potential_energy() 
    
    # This is the optimization routine:
    sopm            =   fmin(optim_cell, 1)[0]
    print 'Static sopm = %.3f' %sopm
    
    
    # BEGIN DYNAMICS
    
    # scales for unit-cell are used to find optimal unit cell in given temperature
    scales      =   [0.98, 1.0, 1.02]
    
    # Store the initial optimal cell
    cell_opm    =   atoms.get_cell().copy()
    
    # Perform dynamics with different cell sizes    
    for i, s in enumerate(scales):
        print s
        p_trajF     =   path_md + 'md_flat_s=%.2f.traj' %s
        p_logF      =   path_md + 'md_flat_s=%.2f.log' %s
        run_dynamics(atoms, opt, params, T, p_trajF, p_logF, 'flat', \
                     s, cell_opm, init_positions, init_cell)
        
    
    
    # Load the trajectories back and get energies, fit parabel and get optimal scale:
    eMean           =   np.zeros(len(scales))
    
    # Load trajectory of each scale:
    for i, s in enumerate(scales):
        p_trajF     =   path_md + 'md_flat_s=%.2f.traj' %s
        
        # Take only thermalized part:
        traj        =   PickleTrajectory(p_trajF)[therm_step:]
        
        # Get the potential energy at each moment:
        eFlatEset   =   np.zeros(len(traj))
        for j, atoms_t in enumerate(traj):
            eFlatEset[j]    =   atoms_t.get_potential_energy()    
        
        # Calculate the average of energy:
        eMean[i]    =   np.mean(eFlatEset)

    # Fit parabel to (scale, eMean) points:        
    
    def x2(sset, *args):
        a, b, c     =   args
        squared     =   np.zeros(len(sset))
        for i, s in enumerate(sset):
            squared[i] = a*s**2 + b*s + c 
        return squared
    
    # Fit the curve, and obtain optimal parameters:    
    c_opm    =   curve_fit(x2, scales, eMean, [1,1,1])[0]
    
    # Optimal scale is then:
    if c_opm[0] > 0:
        # The minima of parabel:
        sopm2    =   -c_opm[1]/(2*c_opm[0])
    else:
        # If no minima -> raise
        print '\nThe appears to be no optimal cell size :P'
        print 'The obtained parabel is upside down -> cant get minima.\n' 
        raise
    
    # Plot this and save figure:
    f, ax       =   plt.subplots(figsize = (4,5))
    dense_scale =   np.linspace(scales[0]*0.97, scales[2]*1.03, 200)
    
    # PLot the thermalized mean energy for each scale:
    ax.scatter(scales, eMean, label='EperA for scales')
    
    # Plot the parabel fit to the thermalized energies
    ax.plot(dense_scale, x2(dense_scale, c_opm[0], c_opm[1], c_opm[2]), label='fit')
    
    # Mark the minima with a dot, optimal scale is then sopm*sopm2 (w.r.t. the initial unit cell)
    ax.scatter(sopm2, x2([sopm2], c_opm[0], c_opm[1], c_opm[2]), label = 'd = %.6f' %(init_cell[0,0]/uzsize[0]*sopm*sopm2))
    
    plt.legend(frameon = False)
    plt.savefig(path_opm + 'scales.svg')
    
    
    # run dynamics with the optimal scale
    print '\nRunning the final dynamics, for optimized unit_cell'
    print 'Optimal scale, in T=%.0fK seems to be s = %.5f,' %(T, sopm*sopm2)
    print 'giving average bond distance d = %.5fAngst. \n' %(init_cell[0,0]/uzsize[0]*sopm*sopm2)
    
    p_trajF     =   path_md + 'md_flat.traj'
    p_logF      =   path_md + 'md_flat.log' 
    run_dynamics(atoms, opt, params, T, p_trajF, p_logF, 'flat', \
                 sopm2, cell_opm, init_positions, init_cell)

    