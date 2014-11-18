'''
Created on 17.11.2014

@author: tohekorh
'''

from ase import *
from hotbit import *
from numpy import *
from ase.lattice.surface import fcc111
from ase import io
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.velocitydistribution import Stationary
from ase.md.langevin import Langevin
from ase.io import PickleTrajectory
from ase import units
from ase.md import MDLogger

T = 100. # temperature
folder = 'data'
nk = 1
r0list = linspace(20,300,50)
fric = 0.002        # equilibration time tau = 10 fs/fric
mdsteps = 1000

d = 2.934
for R in r0list:
    print 'radius',R
    name=(folder,'%.1f' %R)

    atoms = fcc111('Au',size=(2,2,1),a=sqrt(2)*d)
    L = atoms.get_cell().diagonal()
    atoms. set_cell(L) 
    
    angle = L[1]/R
    atoms.rotate('y',pi/2)
    atoms.translate( (-atoms[0].x,0,0) )
    atoms.translate((R,0,0))
    atoms = Atoms(atoms=atoms,container='Wedge')
    atoms.set_container(angle=angle,height=L[0],physical=False,pbcz=True)
    
    calc = Hotbit(SCC=False,txt='%s/optimization_%s.cal' %name,kpts=(nk,1,nk),physical_k=False)
    atoms.set_calculator(calc)

    traj    =   PickleTrajectory( 'md_R%.3f-T%.0f.traj' %(R,T),'w',atoms)   
    dyn     =   Langevin(atoms,dt*units.fs,units.kB*T,fric)
    dyn.attach(MDLogger(dyn, atoms, path + 'md_R%.3f-T%.0f.log' %(R,T), header=False, stress=False,
           peratom=True, mode="w"), interval = 1)
    dyn.run(mdsteps)
    
   
    cohesion = mean([atoms.get_potential_energy() for a in traj]) /len(atoms)
    
    savez('%s/data_%s.npz' %name,cohesion=cohesion)
    io.write('%s/extended_%s.traj' %name,atoms2)