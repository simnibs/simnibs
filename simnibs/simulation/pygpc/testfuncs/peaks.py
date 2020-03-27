# -*- coding: utf-8 -*-
import numpy as np

def peaks(x):
    """ 2-dimensional peaks function.
    
        Input:  
                x       ... input data [N_input x 2]
        Output: 
                res     ... result [N_input x 1]
    """
    y = np.array([(3.0*(1-x[:,0])**2.*np.exp(-(x[:,0]**2) - (x[:,1]+1)**2) \
        - 10.0*(x[:,0]/5.0 - x[:,0]**3 - x[:,1]**5)*np.exp(-x[:,0]**2-x[:,1]**2)\
        - 1.0/3*np.exp(-(x[:,0]+1)**2 - x[:,1]**2))]).transpose()
    
    return y