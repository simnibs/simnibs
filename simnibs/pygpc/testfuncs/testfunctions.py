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