# -*- coding: utf-8 -*-
import numpy as np

def welch_1992(x):
    """ 20-dimensional test function of WELCH (1992)
        
        For input variable screening purposes, it can be found that some input 
        variables of this function have a very high effect on the output, 
        compared to other input variables. As Welch et al. (1992) point out, 
        interactions and nonlinear effects make this function challenging. 
        
        Welch, W. J., Buck, R. J., Sacks, J., Wynn, H. P., Mitchell, T. J., & 
        Morris, M. D. (1992). Screening, predicting, and computer experiments. 
        Technometrics, 34(1), 15-25.
        
        y = welch_1992(x)
        
        Input:  
                x ... input data [N_input x 20] 
                      xi ~ U(-0.5, 0.5), for all i = 1, …, 20.
                
        Output: 
                y ... result [N_input x 1]
    """
    y = np.array([(5.0*x[:,11]/(1+x[:,0]) + 5*(x[:,3]-x[:,19])**2 + x[:,4] + 40*x[:,18]**3 \
        + 5*x[:,18] + 0.05*x[:,1] + 0.08*x[:,2] - 0.03*x[:,5] + 0.03*x[:,6] \
        - 0.09*x[:,8] - 0.01*x[:,9] - 0.07*x[:,10] + 0.25*x[:,12]**2 - 0.04*x[:,13] \
        + 0.06*x[:,14] - 0.01*x[:,16] - 0.03*x[:,17])]).transpose()
    
    return y    