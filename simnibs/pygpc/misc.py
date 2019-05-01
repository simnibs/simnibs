# -*- coding: utf-8 -*-
import numpy as np
import scipy.special
'''
 pygpc software framework for uncertainty and sensitivity
 analysis of complex systems. See also:
    # https://github.com/konstantinweise/pygpc
    #
    # Copyright (C) 2017-2019 the original author (Konstantin Weise),
    # the Max-Planck-Institute for Human Cognitive Brain Sciences ("MPI CBS")
    # and contributors
    #
    # This program is free software: you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation, either version 3 of the License, or
    # (at your option) any later version.
    #
    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.
    #
    # You should have received a copy of the GNU General Public License
    # along with this program.  If not, see <https://www.gnu.org/licenses/>
'''

"""
Generate a cartesian product of input arrays.

Parameters
----------
arrays :    list of array-like 1-D arrays to form the cartesian product of
out :       ndarray

Returns
-------
out : ndarray
      2-D array of shape (M, len(arrays)) containing all combinations

Examples
 --------
>>> combvec(([1, 2, 3], [4, 5], [6, 7]))
array([[1, 4, 6],
       [1, 4, 7],
       [1, 5, 6],
       [1, 5, 7],
       [2, 4, 6],
       [2, 4, 7],
       [2, 5, 6],
       [2, 5, 7],
       [3, 4, 6],
       [3, 4, 7],
       [3, 5, 6],
       [3, 5, 7]])

"""
def combvec(arrays, out=None):

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        combvec(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out
    
def multi_delete(list_, idx):
    indexes = sorted(idx, reverse=True)
    for index in indexes:
        del list_[index]
    return list_    
    
def unique_rows(data):
    uniq, idx = np.unique(data.view(data.dtype.descr * data.shape[1]),return_index=True)
    return uniq[np.argsort(idx)].view(data.dtype).reshape(-1, data.shape[1])
   
def vchoosek(v,k):
    # v is numpy array [1 x nv]
    # k is scalar
    nv = v.shape[1]
    if nv <= k:         
        if nv == k:
            return v
        else:
            return []
    
    d  = nv - k
    ny = d + 1;
    
    for i in range(2,k+1):
        ny = ny + (1.0 * ny * d) / i;
    
    y = np.zeros([int(ny), int(k)])
    
    Index = np.linspace(1,k,k).astype(int)
    Limit = np.append(np.linspace(nv-k+1, nv-1, k-1), 0)
    a = int(1)
    
    while 1:
        b = int(a + nv - Index[k-1])                    # Write index for last column
        for i in range(1,k):                            # Write the left K-1 columns
            y[(a-1):(b), i-1] = v[0,Index[i-1]-1]
         
        y[(a-1):(b), k-1] = v[0,(Index[k-1]-1):(nv)]    # Write the K.th column
        a = b + 1;                                      # Move the write pointer
         
        newLoop = np.sum(Index < Limit)
        if newLoop == 0:             # All columns are filled:
           break                        # Ready!
        Index[(newLoop-1):(k)] = Index[newLoop-1] + np.linspace(1,k-newLoop+1,k-newLoop+1)
    
    return y 
    
def allVL1leq(n, L1):
    res = []
    for i_L1 in range(L1+1):
        # Chose (n-1) the splitting points of the array [0:(n+L1)]
        s = vchoosek(np.linspace(1,n+i_L1-1,n+i_L1-1)*np.ones([1,n+i_L1-1]),n-1) # 1:n+L1-1
        
        m = s.shape[0]

        s1 = np.zeros([m,1])
        s2 = (n+i_L1)+s1

        v = np.diff(np.hstack([s1, s, s2])) 
        v = v-1

        if i_L1 == 0:         
            res = v
        else:
            res = np.vstack([res,v])
        
    return res.astype(int)

def NRMSD(x,ref):
    """ Determine normalized root mean square deviation between input data and
        reference data in [%]

        Input:
                x   ... ND np.array of test data         [ (x), y0, y1, y2 ... ]
                        
                ref ... 1D/2D np.array of reference data [ (xref), yref ]
                        When ref is 1D, all sizes have to match  
                
        Output:
                NRMSD ... normalized root mean square deviation
    """
    N_points = x.shape[0]
    
    # determine number of input arrays
    if ref.shape[1] == 2:
        N_data = x.shape[1]-1
    else:
        N_data = x.shape[1]
    
    # interpolate x on ref data if necessary    
    if ref.shape[1] == 1:
        data = x
        data_ref = ref
    else:
        # crop reference if it is longer than the axis of the data
        xref = ref[ (ref[:,0]>=min(x[:,0])) & (ref[:,0]<=max(x[:,0]))  ,0]
        data_ref = ref[ (ref[:,0]>=min(x[:,0])) & (ref[:,0]<=max(x[:,0])) ,1]
        
        data = np.zeros([len(xref),N_data])
        for i_data in range(N_data):
            data[:,i_data] = np.interp(xref,x[:,0],x[:,i_data+1])
    
    
    if (max(data_ref) - min(data_ref)) == 0: 
        delta = max(data_ref)
    else:
        delta = max(data_ref) - min(data_ref)
    
    # determine NRMSD
    NRMSD = np.zeros(N_data)    
    for i_data in range(N_data):
        NRMSD[i_data] = 100 * np.sqrt( 1.0/N_points * np.sum((data[:,i_data] - data_ref)**2,axis=0)) / delta                         
    
    return NRMSD            
    
def fit_betapdf(data, BETATOL=0, PUNI=0, PLOT=0, VISI=1, filename=[], xlabel="$x$", ylabel="$p(x)$"):
    """ fit_betapdf fits data to a beta distribution in the inverall [a, b] 
    
     fit_betapdf(data, BETATOL=0, PUNI=0, PLOT=1, VISI=1, filename=[], xlabel="$x$", ylabel="$p(x)$"):
    
    Input:
        data    ... data to fit pdf (np.array)
        BETATOL ... tolerance interval to calculate the bounds of beta 
                    distribution from observed data, e.g. 0.2 (+-20% tolerance)
        PUNI    ... [0...1] uniform distribution interval defined as fraction of
                    beta distribution interval, e.g. 0.90 (90%) (default = 0)
        PLOT    ... [1/0] plot fitted pdfs (default = 0)
        VISI    ... [1/0] DefaultFigureVisible on/off (default = 1)
        filename ... filename of plotted data
    
    Output:
        beta_paras ... [p, q, a, b] (2 shape parameters and limits)
        moments    ... [data_mean, data_std, beta_mean, beta_std]
        p_value    ... p-value of the Kolmogorov Smirnov test
        uni_paras  ... [a, b] (limits)
    """
    import matplotlib.pyplot as plt
    data_mean = np.mean(data)
    data_std = np.std(data)
    
    # fit beta distribution to data
    if BETATOL > 0:
        # use user BETATOL of to set limits of distribution
        data_range = data.max()-data.min()
        a_beta = data.min()-BETATOL*data_range
        b_beta = data.max()+BETATOL*data_range
        p_beta, q_beta, a_beta, ab_beta = scipy.stats.beta.fit(data, floc=a_beta, fscale=b_beta-a_beta)
    else:
        # let scipy.stats.beta.fit determine the limits
        p_beta, q_beta, a_beta, ab_beta = scipy.stats.beta.fit(data)
        b_beta = a_beta + ab_beta
    
    beta_mean, beta_var = scipy.stats.beta.stats(p_beta, q_beta, loc=a_beta, scale=(b_beta-a_beta), moments='mv')
    beta_std = np.sqrt(beta_var)
    
    moments = np.array([data_mean, data_std, beta_mean, beta_std])
    
    # determine limits of uniform distribution [a_uni, b_uni] covering the 
    # interval PUNI of the beta distribution
    if PUNI>0:
        a_uni = scipy.stats.beta.ppf((1-PUNI)/2, p_beta, q_beta, loc=a_beta, scale=b_beta-a_beta)
        b_uni = scipy.stats.beta.ppf((1+PUNI)/2, p_beta, q_beta, loc=a_beta, scale=b_beta-a_beta)
        
    # determine kernel density estimates using Gaussian kernel
    kde_x = np.zeros(100)
    kde_y = np.zeros(100)
        
    kde   = scipy.stats.gaussian_kde(data, bw_method=0.05/data.std(ddof=1))
    kde_x = np.linspace(a_beta, b_beta, 100)
    kde_y = kde(kde_x)
    
    # perform Kolmogorov Smirnov test
    D, p_value = scipy.stats.kstest(data, "beta", [p_beta, q_beta, a_beta, ab_beta])
    
    # plot results    
    if PLOT:
        
        if VISI==0:
            plt.ioff()
        else:
            plt.ion()
            
        plt.figure(1) # , figsize=[6, 4]
        plt.clf()
        plt.rc('text', usetex=True)
        plt.rc('font', size=18)
        ax = plt.gca()
        legendtext = [r"e-pdf", r"$\beta$-pdf"]
        
        # plot histogram of data
        n, bins, patches = plt.hist(data, bins=10, normed=1, color=[1,1,0.6], alpha=0.5)
        
        # plot beta pdf (kernel density estimate)
        plt.plot(kde_x, kde_y, 'r--', linewidth=2)  
        
        # plot beta pdf (fitted)
        beta_x = np.linspace(a_beta, b_beta, 100)
        beta_y = scipy.stats.beta.pdf(beta_x,p_beta,q_beta,loc=a_beta,scale=b_beta-a_beta)
        
        plt.plot(beta_x, beta_y, linewidth=2, color = [0, 0, 1])
        
        # plot beta pdf (kernel density estimate)
        uni_y = 0
        if PUNI > 0:
            uni_x = np.hstack([a_beta, a_uni-1E-6*(b_uni-a_uni), np.linspace(a_uni,b_uni,100),b_uni+1E-6*(b_uni-a_uni), b_beta])
            uni_y = np.hstack([0, 0, 1.0/(b_uni-a_uni)*np.ones(100), 0, 0])
            plt.plot(uni_x, uni_y, linewidth=2, color = 'k') 
            legendtext.append("u-pdf")               
        
        # configure plot
        plt.legend(legendtext, fontsize=18, loc = "upper left")         
        plt.grid(True) 
        plt.xlabel(xlabel,fontsize=22)
        plt.ylabel(ylabel,fontsize=22)
        ax.set_xlim(a_beta-0.05*(b_beta-a_beta),b_beta+0.05*(b_beta-a_beta))
        ax.set_ylim(0,1.15*max([max(data),max(beta_y),max(uni_y)]))
        
        if VISI>0:
            plt.show()
            
        # save plot
        if filename:
            plt.savefig(filename + ".pdf", format='pdf', bbox_inches='tight', pad_inches=0.01*4)
            plt.savefig(filename + ".png", format='png', bbox_inches='tight', pad_inches=0.01*4, dpi=600)
    
    # return results
    beta_paras = np.array([p_beta, q_beta, a_beta, b_beta])
    if PUNI > 0:
        uni_paras  = np.array([a_uni, b_uni])
        return beta_paras, moments, p_value, uni_paras
    else:
        return beta_paras, moments, p_value     
            
def mutcoh(A):
    """ mutcoh calculates the mutual coherence of a matrix A, i.e. the cosine
        of the smallest angle between two columns
      
        mutual_coherence = mutcoh(A)
 
        Input:
            A ... m x n matrix
 
        Output:
            mutual_coherence """

    T = np.dot(A.conj().T,A)
    s = np.sqrt(np.diag(T))    
    S = np.diag(s)
    mutual_coherence = np.max(1.0*(T-S)/np.outer(s,s))
    
    return mutual_coherence
            
        
            
    
    
    
