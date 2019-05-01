# -*- coding: utf-8 -*-
import numpy as np

def ishigami(x,a,b):
    """ 3-dimensional test function of Ishigami.
        The Ishigami function of Ishigami & Homma (1990) is used as an example
        for uncertainty and sensitivity analysis methods, because it exhibits
        strong nonlinearity and nonmonotonicity. It also has a peculiar
        dependence on x3, as described by Sobol' & Levitan (1999).
        
        Ishigami, T., & Homma, T. (1990, December). An importance quantification
        technique in uncertainty analysis for computer models. In Uncertainty
        Modeling and Analysis, 1990. Proceedings., First International Symposium
        on (pp. 398-403). IEEE.
        
        Sobol', I. M., & Levitan, Y. L. (1999). On the use of variance reducing
        multipliers in Monte Carlo computations of a global sensitivity index.
        Computer Physics Communications, 117(1), 52-61.
        
        f(x) = sin(x1) + a*sin(x2)^2 + b*x3^4*sin(x1)        
        
        y = ishigami(x,a,b)
        
        Input:  
                x   ... input data [N_input x 3]
                        xi ~ Uniform[-π, π], for all i = 1, 2, 3
                a,b ... shape parameter
                
        Output: 
                y   ... result [N_input x 1]
    """
    y = np.array([(np.sin(x[:,0])+a*np.sin(x[:,1])**2+b*x[:,2]**4*np.sin(x[:,0]))]).transpose()
    
    return y