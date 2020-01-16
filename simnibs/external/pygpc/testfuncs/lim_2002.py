# -*- coding: utf-8 -*-
import numpy as np

def lim_2002(x):
    """ 2-dimensional test function of Lim et al.
        This function is a polynomial in two dimensions, with terms up to degree
        5. It is nonlinear, and it is smooth despite being complex, which is
        common for computer experiment functions (Lim et al., 2002). 
        
        Lim, Y. B., Sacks, J., Studden, W. J., & Welch, W. J. (2002). Design
        and analysis of computer experiments when the output is highly correlated
        over the input space. Canadian Journal of Statistics, 30(1), 109-126.
        
        f(x) = 9 + 5/2*x1 - 35/2*x2 + 5/2*x1*x2 + 19*x2^2 - 15/2*x1^3 
               - 5/2*x1*x2^2 - 11/2*x2^4 + x1^3*x2^2         
        
        y = lim_2002(x)
        
        Input:  
                x ... input data [N_input x 2]
                      xi ∈ [0, 1], for all i = 1, 2
                
        Output: 
                y ... result [N_input x 1]
    """
    
    y = np.array([(9 + 5.0/2*x[:,0] - 35.0/2*x[:,1] + 5.0/2*x[:,0]*x[:,1] 
        + 19*x[:,1]**2 - 15.0/2*x[:,0]**3 - 5.0/2*x[:,0]*x[:,1]**2 - 11.0/2*x[:,1]**4 
        + x[:,0]**3*x[:,1]**2)]).transpose()
    
    return y