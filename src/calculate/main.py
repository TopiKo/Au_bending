'''
Created on 17.11.2014

@author: tohekorh
'''
import numpy as np
from sphereDyn import sphereDynamics
from cylinderDyn import cylinderDynamics
from flatDyn import flatDynamics
#from scipy.optimize import curve_fit
from aid import checkIfExists, get_Rset

# Root path for the defferent files:
path        =   '/space/tohekorh/Au_bend/files/'

# The temperatures used in the simulation:
TList       =   [12]

# The unit cell angles, preferring 2*pi/integer
angles      =   [np.pi/2, np.pi/6, np.pi/18, np.pi/36]

# Langevin, friction and time step, and simul time = (1+m)*relax_time 
fric        =   0.002
dt          =   5.0
m           =   3

# Unit cell size (8,4,1) = square 64 atoms, bond distance quess = d.
uz_size     =   (8,4,1)
d           =   2.934


params              =   {}

params['d']         =   d           # lattice constant for gold
params['nk']        =   3           # max number of k-points
params['fmax']      =   1E-2        # Max force for relaxation
params['fric']      =   fric        # equilibration time tau = 10 fs/fric (= 1/M*mdsteps*dt) 
params['dt']        =   dt          # Timestep for moldy
tau                 =   10./fric    # This is the thermalization time

simul_time          =   tau + m*tau # The simulation time is thermalization time plus some time after thermalization 

params['mdsteps']   =   9 #int(simul_time/dt)          # Number of md steps
params['thermN']    =   params['mdsteps']/( m + 1 ) # From this step onward system is supposed to be thermalized
params['path']      =   path                        
params['uz_size']   =   uz_size


# Loop over temperatures
for iT, T in enumerate(TList):
    print 'temperature = %i' %(int(T))
    
    # Check if the calculation is performed already:
    if not checkIfExists(T, params['mdsteps'], 'flat', path):    
        print 'Flat Dynamics'
        flatDynamics(T, params)
        
    # Get proper radiuces according to the optimized unit-cell 
    # size in temperature T from flat simulation:
    r0list_c    =   get_Rset(angles, T, params)


    # CYLINDER DYNAMICS
    for ir, R in enumerate(r0list_c):

        if not checkIfExists(T, params['mdsteps'], 'cylinder', path, R): 
            print '\nCylinder Dyn'
            cylinderDynamics(R, T, params)
            
    
    # Forget this, for now...
    # SPHERE DYNAMICS
    #for ir, R in enumerate(r0list_s):

    #    if not checkIfExists(T, params['mdsteps'], 'sphere', path, R): 
    #        print 'Sphere Dyn'
    #        sphereDynamics(R, T, params)

    
