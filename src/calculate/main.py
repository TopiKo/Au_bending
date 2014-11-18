'''
Created on 17.11.2014

@author: tohekorh
'''
import numpy as np
from sphereDyn import sphereDynamics
from cylinderDyn import cylinderDynamics
from flatDyn import flatDynamics
#from scipy.optimize import curve_fit



r0list  =   np.linspace(20,100,8)
rlist_L =   np.linspace(r0list[0] - 3, r0list[-1] + 10, 200)

path    =   '/space/tohekorh/Au_bend/files/'

TList   =   np.arange(500, 2101, 200) 
fric    =   0.002
dt      =   2.0

params  =   {}

params['d']     =   2.934
params['nk']    =   20
params['fmax']  =   1 #1E-2
params['fric']  =   fric        # equilibration time tau = 10 fs/fric (= mdsteps*dt) 
params['dt']    =   dt
tau             =   40./fric  

params['mdsteps']   =   1 #int(tau/dt)
params['path']      =   path


print params['mdsteps']


for iT, T in enumerate(TList):
    print 'temperature = %i' %(int(T))
    flatDynamics(T, params)
    
    for ir, R in enumerate(r0list):
        cylinderDynamics(R, T, params)
        sphereDynamics(R, T, params)
        
