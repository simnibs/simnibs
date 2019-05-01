# -*- coding: utf-8 -*-
import numpy as np

def sphere0_fun(x,a,b):
    """ N-dimensional sphere function with zero mean.
        
        y = sphere0(x,a,b)

        Input:  
                x   ... input data [N_input x N_dims]
                a,b ... lower and upper bound of all input vars
        Output: 
                y   ... result [N_input x 1]
    """

    try:
        N = x.shape[1]
    except IndexError:
        N=1
        x=np.array([x])

    # zero mean   
    c2 = (1.0*N*(b**3-a**3))/(3*(b-a))
    
    # sphere function
    y = np.array([(np.sum(np.square(x),axis=1)-c2)]).transpose()
    
    return y  