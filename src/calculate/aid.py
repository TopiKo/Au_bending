'''
Created on 18.11.2014

@author: tohekorh
'''
import os

def checkAndCreateFolder(path):
    
    if not os.path.exists(path):
        os.makedirs(path)