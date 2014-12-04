'''
Created on 27.11.2014

@author: tohekorh
'''
import numpy as np
from ase import *
from hotbit import *
from ase.visualize import view
import matplotlib.pyplot as plt 
from ase.optimize import BFGS
from aid import get_square_uCell

path        =   '/space/tohekorh/Au_bend/files/'

def test(size, R, nk):

    
    atoms_flat  =   get_square_uCell(size)
    
    view(atoms_flat)
    
    # CALCULATOR FLAT
    calc_f      =   Hotbit(SCC=False, kpts=(nk,nk,1),  \
                           txt= path + 'test_consistency/optimization_flat.cal')
    atoms_flat.set_calculator(calc_f)
    
    opt_f       =   BFGS(atoms_flat)
    opt_f.run(fmax = 0.05)
    e_flat      =   atoms_flat.get_potential_energy()
    

    
    
    atoms_c     =   atoms_flat.copy()
    
    L           =   atoms_c.get_cell().diagonal()
     
    atoms_c.set_cell(L) 
    angle       =   L[1]/R
    atoms_c.rotate('y', np.pi/2)
    atoms_c.translate((-atoms_c[0].x, 0, 0) )
    
    for a in atoms_c:
        r0      =   a.position
        phi     =   r0[1]/L[1]*angle
        a.position[0]   =   R*np.cos(phi)
        a.position[1]   =   R*np.sin(phi)
    
    
    atoms_c       =   Atoms(atoms = atoms_c, container = 'Wedge')
    atoms_c.set_container(angle = angle, height = L[0], physical = False, pbcz = True)
    
    if R < 100:
        view(atoms_c.extended_copy((8,1,3)))
    
    # CALCULATOR Cyl
    calc_c      =   Hotbit(SCC=False, kpts=(nk,1, nk), physical_k = False, \
                           txt= path + 'test_consistency/optimization_cyl.cal')
    atoms_c.set_calculator(calc_c)
    
    opt_c       =   BFGS(atoms_c)
    opt_c.run(fmax = 0.05)
    
    
    e_cyl       =   atoms_c.get_potential_energy()

    print 'R = %.2f' %R
    print 'energy flat     = %.6f' %e_flat
    print 'energy cylinder = %.6f' %e_cyl
    print 'energy dif (e_cylinder - eflat)/nAtoms  = %.6f \n' %(-(e_flat - e_cyl)/len(atoms_flat))
    
    return e_flat, e_cyl, len(atoms_flat)
    
    
Rs      =   [60, 240, 5000, 1e4]    
uzS     =   [(1,1,1), (2,2,1), (3,3,1)]
nks     =   [6, 3, 2]
data    =   np.zeros((len(Rs), 3))
f, ax   =   plt.subplots(figsize = (4,5))

for k, uz in enumerate(uzS):
    print 'uzsize = ', uz
    for i, r in enumerate(Rs):
        ef, ec, nAtoms  =   test(uz, r, nks[k])
        data[i]         =   [r, ef/nAtoms, ec/nAtoms]
    
    
    ax.plot(data[:,0], data[:,1], '-o', label='eF, Natoms_uz = %.i' %nAtoms)
    ax.plot(data[:,0], data[:,2], '-o', label='eC, Natoms_uz = %.i' %nAtoms)
    ax.plot(data[:,0], data[:,2] - data[:,1], '-o', label='EC - eF, Natoms_uz = %.i' %nAtoms)
    
ax.set_xlabel('radius')
ax.set_ylabel('e deviation per atom')
plt.legend()
plt.show()

#f, ax   =   plt.subplots(figsize = (4,5))
#ax.plot(dev[:,0], dev[:,1], '-o', label='dev')
#ax.set_xlabel('uzell size')
#ax.set_ylabel('e deviation per atom')

plt.show()