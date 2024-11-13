# -*- coding: utf-8 -*-
import numpy as np
   
def oakley_ohagan_2004(x):
    """ 15-dimensional test function of OAKLEY & O'HAGAN (2004)
        
        This function's a-coefficients are chosen so that 5 of the input
        variables contribute significantly to the output variance, 5 have a 
        much smaller effect, and the remaining 5 have almost no effect on the 
        output variance. 
        
        Oakley, J. E., & O'Hagan, A. (2004). Probabilistic sensitivity analysis 
        of complex models: a Bayesian approach. Journal of the Royal Statistical
        Society: Series B (Statistical Methodology), 66(3), 751-769.
        
        y = oakley_ohagan_2004(x)
        
        Input:  
                x ... input data [N_input x 15] 
                      xi ~ N(μ=0, σ=1), for all i = 1, …, 15.
                
        Output: 
                y ... result [N_input x 1]
    """
    
    # load coefficients
    M = np.loadtxt('misc/oakley_ohagan_2004_M.txt')
    a1 = np.loadtxt('misc/oakley_ohagan_2004_a1.txt')
    a2 = np.loadtxt('misc/oakley_ohagan_2004_a2.txt')
    a3 = np.loadtxt('misc/oakley_ohagan_2004_a3.txt')
    
    # function
    y = np.array([(np.dot(x,a1) + np.dot(np.sin(x),a2) + np.dot(np.cos(x),a3) \
        + np.sum(np.multiply(np.dot(x,M),x),axis=1))]).transpose()
    
    return y