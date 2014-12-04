'''
Created on 17.11.2014

@author: tohekorh
'''
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt 
from ase.io import PickleTrajectory
from aid import checkIfExists, trim, get_Rset,  checkAndCreateFolder, get_AvRadius, get_Area
from ase.visualize import view

# This plotting module

# Path is the root for the files as defined in main.py
path        =   '/space/tohekorh/Au_bend/files/'

# The angles should match the ones one has calculated:
angles      =   [np.pi/2, np.pi/6, np.pi/18, np.pi/36] 

# The amount of moldy steps in the trajevctories:
nm          =   4000 #2000 #number of moldy steps

# List of temperatures:
TList       =   [1,20,100] #[2, 15, 55] #np.arange(500, 2101, 200)

# From this step onwart the system is supposed to be thermalized
therm       =   int(nm/4.)

# If temp = False do the static case = read 0 temperature optimized structures work with them.
# If temp = True use the whole thermalized trajectories and treat energies as averages of the 'structures'.
temp        =   False


params              =   {}
params['path']      =   path 
params['mdsteps']   =   nm 






def CoverR2(Rset, *args):
    k      =   args[0]
    return .5*k/Rset**2

c_coll          =   np.zeros((len(TList), 2))

for iT, T in enumerate(TList):
    cont            =   False
    print 'temperature = %i' %(int(T))
    
    if checkIfExists(T, nm, 'flat', path):
        traj_flat   =   PickleTrajectory(path + 'flat/md_data/T=%.0f/md_flat.traj' %T)
        log_flat    =   np.loadtxt(path + 'flat/md_data/T=%.0f/md_flat.log' %T)
        
        
        if temp:
            # Take only configs from thermalization and onward
            traj_flat   =   traj_flat[therm:]
            log_flat    =   log_flat[therm:]
            
            # Get area of the atoms (=area of the unit cell)
            A           =   get_Area(traj_flat)
            natoms      =   len(traj_flat[0])
            
            # log_flat[:,1] is the toal energy of the system:
            eflat_table =   log_flat[:,1]*natoms/A
            
            eflatPerArea       =   np.mean(eflat_table)
            
        else:
            # Consider only the otimized structure:
            traj_flat   =   PickleTrajectory(path + 'flat/opm/T=%.0f/optimization_flat.traj' %T)
            
            # Last one is the fully optimized:
            atoms_flat  =   traj_flat[-1]
            
            A           =   get_Area([atoms_flat])
            natoms      =   len(atoms_flat)
            eflatPerArea       =   atoms_flat.get_potential_energy()/A
        
        cont        =   True
    
    else:
        print 'ei oo flat'
    
            
    
    if cont:
        # Get the set of radiuces that should have been calculated wwith the give angles:
        r0list_c    =   get_Rset(angles, T, params)
        
        # Dense list for plottin:
        rlist_L     =   np.linspace(r0list_c[0]*0.99, r0list_c[-1]*1.01, 200)
        
        # Table to save stuff:
        bEnergyPerArea_c = np.zeros((len(r0list_c),3))

        for ir, R in enumerate(r0list_c):

            if checkIfExists(T, nm, 'cylinder', path, R):
                traj_c  =   PickleTrajectory(path + 'cylinder/md_data/T=%.0f/md_R%.3f.traj' %(T, R))
                log_c   =   np.loadtxt(path + 'cylinder/md_data/T=%.0f/md_R%.3f.log' %(T, R))

                if temp:
                    # Take only configs from thermalization and onward
                    traj_c      =   traj_c[therm:]
                    log_c       =   log_c[therm:]

                    Rav, phi, h =   get_AvRadius(traj_c)
                    natoms      =   len(traj_c[0])
                    A_w         =   Rav*phi*h
                    ec_set      =   natoms*log_c[:,1]/A_w
                
                    
                    # check if the energies are normally distributed
                    e_division      =   np.linspace(min(ec_set), max(ec_set), 10)
                    n_division      =   np.zeros(len(e_division))
                                
                    for i in range(len(e_division) - 1):
                        emin = e_division[i]
                        emax = e_division[i+1] 
                        
                        for e in ec_set:
                            if emin < e < emax:
                                n_division[i] += 1
                           
                    f, ax           =   plt.subplots(figsize = (4,4))
                    ax.plot(e_division, n_division)
                    plt.show()
                    # End check
                    
                    ec          =   np.mean(ec_set) - eflatPerArea
                    std_dev     =   np.std(ec_set)
                    mean_err    =   std_dev/np.sqrt(len(ec_set))
                    bEnergyPerArea_c[ir]  =   [Rav, ec, mean_err]   

                else:
                    traj_c      =   PickleTrajectory(path + 'cylinder/opm/T=%.0f/optimization_%.1f.traj' %(T, R))
                    traj_c2     =   PickleTrajectory(path + 'cylinder/md_data/T=%.0f/md_R%.3f.traj' %(T, R))

                    atoms_c     =   traj_c[-1]
                    
                    Rav, phi, h =   get_AvRadius([atoms_c])
                    natoms      =   len(atoms_c)
                    A_w         =   Rav*phi*h
                    ec_set      =   atoms_c.get_potential_energy()/A_w
                    
                    ec          =   ec_set - eflatPerArea
                    
                    bEnergyPerArea_c[ir]  =   [Rav, ec, 0]
                    
                    # Cheks:
                    #print  Rav, R, phi, h 
                    print  R, atoms_c.get_potential_energy(), traj_c2[0].get_potential_energy() #atoms_flat.get_potential_energy()
                    #print  traj_c[0].get_potential_energy() - len(traj_c[0])*log_c[1]
                    #print  traj_flat[0].get_potential_energy() - len(traj_c[0])*log_flat[1]
                    #view(traj_c[0])
                    #print  ec_set, eflatPerArea
                    #print  ec 
                    #print
               
            else:
                print 'ei oo sylinteria' 
                 
        # Remove all zeros from table (0 = not yet calculated)
        bEnergyPerArea_c = trim(bEnergyPerArea_c)
        
        #print bEnergyPerArea_c
        
        # Fit y = .5*kappa/R**2 
        c_opm, c_cov  =   curve_fit(CoverR2, bEnergyPerArea_c[:,0], bEnergyPerArea_c[:,1], p0 = 1) #C_cov describes the error   
        
        print c_opm
        c_coll[iT]    =   [T, c_opm[0]]
        
        # PLOT
        f, ax1   =   plt.subplots(figsize = (4,5))
        
        plt.rc('text', usetex = True)
        plt.rc('font', family = 'serif')
          
        ax1.errorbar(bEnergyPerArea_c[:,0], bEnergyPerArea_c[:,1], yerr = bEnergyPerArea_c[:,2], \
                    fmt = '--', color = 'black', label='data')
        fit_c   =    CoverR2(rlist_L, c_opm[0]) #, c_opm[1])
        ax1.plot(rlist_L, fit_c, label = 'curve fit, Kappa = %.3f' %(c_opm[0])) 
        
        ax1.set_xlabel(r'Radius [$\AA$]') 
        ax1.legend(frameon = False)
        ax1.set_title('(E(R) - eflatPerArea)/area, T=%i' %(int(T)))
        
        path_pic    =   path + 'cylinder/pictures/' 
        checkAndCreateFolder(path_pic)
        plt.savefig(path_pic + 'cylinder_T=%i.svg' %(int(T)))
        plt.show()
        header = '#radius Angst, Cohesion [eV], mean error Cohesion'
        #np.savetxt(path + 'sphere/cohesion_T=%i.data' %(int(T)), cohesion_s, header = header)

header = '#T[K], Cs_opm[eV], Cs_opmError, c_opm[eV], Cc_opmError'
#np.savetxt(path + 'sphere/C_sphere.data', Cs, header = header)


f, ax   =   plt.subplots(figsize = (8,5))

kappas  =   c_coll[:,1]
    
ax.plot(c_coll[:,0], kappas, '-o', label='kappa')
ax.set_ylabel('Bend modulus Kappa')
ax.set_xlabel('Temperature')
plt.show()

print c_coll
