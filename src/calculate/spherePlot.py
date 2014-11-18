'''
Created on 17.11.2014

@author: tohekorh
'''
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt 
from ase.io import PickleTrajectory


r0list  =   np.linspace(20,100,8)
rlist_L =   np.linspace(r0list[0] - 3, r0list[-1] + 10, 200)

path    =   '/space/tohekorh/Au_bend/files/'

TList   =   np.arange(500, 2101, 200)



def CoverR2(Rset, *args):
    C = args[0]
    return C/Rset**2



Cs  =   np.zeros((len(TList), 3))

for iT, T in enumerate(TList):
    cont            =   False
    print 'temperature = %i' %(int(T))
    try:
        traj_flat   =   PickleTrajectory(path + 'flat/md_data/T=%.0f/md_flat.traj' %T)
        eflat       =   np.mean([atoms.get_potential_energy() / len(atoms) for atoms in traj_flat])
        cont        =   True
    except IOError:
        print 'eiooo'
            
    
    if cont:
        cohesion_s = np.zeros((len(r0list),3))
        cohesion_c = np.zeros((len(r0list),3))

        for ir, R in enumerate(r0list):
            
            try:
                traj_c          =   PickleTrajectory(path + 'cylinder/md_data/T=.0f/md_R%.3f.traj' %(T, R))
                ec_set          =   [atoms.get_potential_energy() / len(atoms) - eflat for atoms in traj_c]
                '''
                # check if the energies are normally distributed
                e_division      =   np.linspace(min(es_set), max(es_set), 200)
                n_division      =   np.zeros(len(e_division))
                            
                for i in range(len(e_division) - 1):
                    emin = e_division[i]
                    emax = e_division[i+1] 
                    
                    for e in es_set:
                        if emin < e < emax:
                            n_division[i] += 1
                       
                f, ax           =   plt.subplots(figsize = (4,4))
                ax.plot(e_division, n_division)
                plt.show()
                '''
                
                ec              =   np.mean(ec_set)
                std_dev         =   np.std(ec_set)
                mean_err        =   std_dev/np.sqrt(len(ec_set))
                cohesion_c[ir]  =   [R, ec, mean_err] 
                
            except IOError:
                print ir
                print 'eiooo'
        
            
            
            try:
                traj_s          =   PickleTrajectory(path + 'sphere/md_data/T=.0f/md_R%.3f.traj' %(T, R))
                es_set          =   [atoms.get_potential_energy() / len(atoms) - eflat for atoms in traj_s]
                '''
                # check if the energies are normally distributed
                e_division      =   np.linspace(min(es_set), max(es_set), 200)
                n_division      =   np.zeros(len(e_division))
                            
                for i in range(len(e_division) - 1):
                    emin = e_division[i]
                    emax = e_division[i+1] 
                    
                    for e in es_set:
                        if emin < e < emax:
                            n_division[i] += 1
                       
                f, ax           =   plt.subplots(figsize = (4,4))
                ax.plot(e_division, n_division)
                plt.show()
                '''
                
                es              =   np.mean(es_set)
                std_dev         =   np.std(es_set)
                mean_err        =   std_dev/np.sqrt(len(es_set))
                cohesion_s[ir]    =   [R, es, mean_err] 
                
            except IOError:
                print ir
                print 'eiooo'
                
        
        # Fit C/R**2 curve to E(R) - E_flat    
        C_opm, C_cov    =   curve_fit(CoverR2, r0list, cohesion_s[:,1], p0 = 1) #C_cov describes the error   
        
        Cs[iT]          =   [T, C_opm, C_cov]
    
        # PLOT
        f, ax           =   plt.subplots(figsize = (4,4))
        legp, leg       =   [], []
        c               =   ['b', 'g', 'r']
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')  
        #ax.scatter(r0list, cohesion[:,1], label = 'Cohesion, T=%i' %(int(T)))
        ax.errorbar(cohesion_s[:,0], cohesion_s[:,1], yerr = cohesion_s[:,2], \
                    fmt = '--', color = 'black', label='Cohesion, T=%i' %(int(T)))
        ax.plot(rlist_L, CoverR2(rlist_L, C_opm), label = 'curve fit, C = %.3f' %C_opm) 
        ax.set_xlabel(r'Radius [$\AA$]') 
        ax.set_ylabel(r'Cohesion = Es - Eflat, [eV]') 
        ax.legend(frameon = False)
        plt.savefig(path + 'sphere/pictures/' + 'sphere_T=%i.svg' %(int(T)))
        plt.show()
        
        header = '#radius Angst, Cohesion [eV], mean error Cohesion'
        np.savetxt(path + 'sphere/cohesion_T=%i.data' %(int(T)), cohesion_s, header = header)

header = '#T[K], C_opm[eV], C_opmError'
np.savetxt(path + 'sphere/C_sphere.data', Cs, header = header)
    
