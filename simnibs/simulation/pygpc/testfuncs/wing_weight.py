# -*- coding: utf-8 -*-
import numpy as np

def wing_weight(x):
    """ 10-dimensional test function which models a light aircraft wing
        
        Forrester, A., Sobester, A., & Keane, A. (2008). Engineering design via 
        surrogate modelling: a practical guide. Wiley.
        
        y  = wing_weight(x)
        
        Input:  
                x ... input data [N_input x 10] 
                      x1(Sw)  ∈ [150, 200]
                      x2(Wfw) ∈ [220, 300]
                      x3(A)   ∈ [6, 10]
                      x4(Λ)   ∈ [-10, 10]
                      x5(q)   ∈ [16, 45]
                      x6(λ)   ∈ [0.5, 1]
                      x7(tc)  ∈ [0.08, 0.18]
                      x8(Nz)  ∈ [2.5, 6]
                      x9(Wdg) ∈ [1700, 2500]
                      x10(Wp) ∈ [0.025, 0.08]
                
        Output: 
                y ... result [N_input x 1]
    """
    y = np.array([( 0.036*x[:,0]**0.758 * x[:,1]**0.0035 \
          * (x[:,2]/np.cos(x[:,3])**2)**0.6 * x[:,4]**0.006 * x[:,5]**0.04 \
          * (100*x[:,6]/np.cos(x[:,3]))**-0.3 * (x[:,7]*x[:,8])**0.49 \
          + x[:,0]*x[:,9])]).transpose()
    
    return y       