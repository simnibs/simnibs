# -*- coding: utf-8 -*-
import numpy as np

def sphere_fun(x):
    """ N-dimensional sphere function with zero mean.
        
        y = sphere0(x)
        
        Input:  
                x ... input data [N_input x N_dims]
        Output: 
                y ... result [N_input x 1]
    """
    
    # sphere function
    y = np.array([(np.sum(np.square(x),axis=1))]).transpose()
    
    return y