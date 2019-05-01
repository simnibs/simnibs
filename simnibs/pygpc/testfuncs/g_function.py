# -*- coding: utf-8 -*-
import numpy as np

def g_function(x,a):
    """ N-dimensional g-function used by Saltelli and Sobol
        this test function is used as an integrand for various numerical 
        estimation methods, including sensitivity analysis methods, because it 
        is fairly complex, and its sensitivity indices can be expressed 
        analytically. The exact value of the integral with this function as an 
        integrand is 1. 
        
        Saltelli, Andrea; Sobol, I. M. (1995): Sensitivity analysis for nonlinear
        mathematical models: numerical experience. In: Mathematical models and
        computer experiment 7 (11), S. 16–28.
        
        y = g_function(x,a)
        
        Input:  
                x ... input data [N_input x N_dims]
                a ... importance factor of dimensions [N_dims]
        Output: 
                y ... result [N_input x 1]
    """
         
    try:
        x.shape[1]
    except IndexError:
        x=np.array([x])

    # g-function
    y = np.array([(np.prod((np.abs(4.0*x-2)+a)/(1.0+a),axis=1))]).transpose()
    
    return res     