'''
Created on 17.11.2014

@author: tohekorh
'''

from ase import *
from hotbit import *
import numpy as np
from ase.optimize import BFGS
#from scipy.optimize import curve_fit
from aid import checkAndCreateFolder, get_square_uCell, scale_cell_flat
from scipy.optimize import fmin
from ase.visualize import view
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

path        =   '/space/tohekorh/Au_bend/files/'

angles      =   [np.pi/2, np.pi/6, np.pi/18, np.pi/36]
fric        =   0.002
dt          =   5.0
uz_size     =   (8,4,1)
d           =   2.934
nk          =   3
fmax        =   1E-2 


    


# START WITH FLAT
atoms       =   get_square_uCell(uz_size)
init_cell   =   atoms.get_cell().copy()

# FOLDERS
path_opm    =   path + 'static/opm/'
checkAndCreateFolder(path_opm)

#view(atoms)

# CALCULATOR
calc        =   Hotbit(SCC=False, kpts=(nk,nk,1),  \
                       txt= path_opm + 'optimization_flat.cal')
atoms.set_calculator(calc)

# RELAX
opt         =   BFGS(atoms, trajectory= path_opm + 'optimization_flat.traj')


def optim_cell(s):
    
    cell        =   scale_cell_flat(init_cell, s)
    atoms.set_cell(cell, scale_atoms = True)
    opt.run(fmax)
    return atoms.get_potential_energy() 

sopm    =   fmin(optim_cell, 1)[0]
print 'Static sopm = %.7f' %sopm

eflatPerArea   =   atoms.get_potential_energy()

L       =   atoms.get_cell().diagonal() 
r0list  =   np.zeros(len(angles))
for i, angle in enumerate(angles):
    r0list[i]   =   L[1]/angle

ec      =   np.zeros(len(r0list))
r_real  =   np.zeros(len(r0list))
# CYLINDER DYNAMICS
for ir, R in enumerate(r0list):
    
    atoms_c     =   atoms.copy()
    angle       =   L[1]/R
    atoms_c.rotate('y', np.pi/2)
    atoms_c.translate((-atoms_c[0].x, 0, 0) )
    
    # Initial map for atoms
    phi_max     =   angle
    for a in atoms_c:
        r0      =   a.position
        phi     =   r0[1]/L[1]*phi_max
        a.position[0]   =   R*np.cos(phi)
        a.position[1]   =   R*np.sin(phi)

    #view(atoms_c)
    
    
    # proper number of kappa points in angle direction
    if (2*np.pi/angle)%1 > 0:
        raise 

    
    atoms_c     =   Atoms(atoms = atoms_c, container = 'Wedge')
    atoms_c.set_container(angle = angle, height = L[0], physical = False, pbcz = True)
    
    # CALCULATOR
    calc        =   Hotbit(SCC=False, kpts=(nk, 1, nk), physical_k = False, \
                    txt= path_opm + 'optimization_cylR=%.3f.cal' %R)
    atoms_c.set_calculator(calc)
    
    # RELAX
    opt         =   BFGS(atoms_c, trajectory= path_opm + 'optimization_cylR=%.3f.traj' %R)
    opt.run(fmax=fmax, steps=1000)            
    
    ec[ir]      =   atoms_c.get_potential_energy()   
    
    for a in atoms_c:
        r0          =   a.position
        r_real[ir] +=   np.sqrt(r0[0]**2 + r0[1]**2)/len(atoms_c)
    
    print R, ec[ir]


eflatPerArea    =   eflatPerArea/(L[0]*L[1])

print eflatPerArea, L[0]*L[1], eflatPerArea/(L[0]*L[1]) 
ecPerArea       =   np.zeros(len(r0list))
print 
for i in range(len(r0list)):
    ecPerArea[i]=   ec[i]/(L[0]*r_real[i]*angles[i])
    #print  r0list[i], r_real[i] 
    print  ec[i], eflatPerArea
    #print  ecPerArea[i], eflatPerArea 
    #print  ecPerArea[i] - eflatPerArea
    print 
         
def CoverR2(rlist, *args):
    k = args[0]
    return .5*k/rlist**2
    

c_opm, c_cov  =   curve_fit(CoverR2, r0list, ecPerArea - eflatPerArea, p0 = 1)
print c_opm
f, ax   =   plt.subplots(figsize = (4,5))
rlist   =   np.linspace(r0list[0], r0list[-1], 100)

ax.scatter(r0list, ecPerArea - eflatPerArea)
ax.plot(rlist, CoverR2(rlist, c_opm))
plt.show()
print r0list
print ecPerArea - eflatPerArea
