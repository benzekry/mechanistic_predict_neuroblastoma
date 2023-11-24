"""
Created on Mon Aug 19 2019

@author: benzekry
"""
import numpy as np
import pdb
def vol2cell(V):
    '''
    Converts volume in mm3 into number of cells assuming 1 mm3 = 1e6 cells
    '''
    return V*1e6
def vol2diam(V):
    '''
    Converts volume into diameter assuming spherical shape
    '''
    D = (6*V/np.pi)**(1/3)
    return D
def diam2cell(D):
    '''
    Converts diameter into number of cells assuming spherical shape and 1 mm3 = 1e6 cells
    '''
    N = 1/6*np.pi*(D**3)*1e6
    return N
def cell2diam(N):
    '''
    Converts number of cells into diameter assuming spherical shape and 1 mm3 = 1e6 cells
    '''    
    D = vol2diam(N*1e-6)
    return D
