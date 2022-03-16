# -*- coding: utf-8 -*-

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

import numpy as np
import scipy.special
import scipy.stats
import pickle         # module for saving and loading object instances
import h5py

from .grid import quadrature_jacobi_1D
from .grid import quadrature_hermite_1D
from .grid import randomgrid

from .misc import unique_rows
from .misc import allVL1leq


def save_gpcobj_pkl(gobj, fname):
    """ saving gpc object including infos about input pdfs, polynomials, grid etc...
    
    save_gpc(gobj, filename)
    
    gobj     ... gpc object to save
    fname ... filename with .pkl extension """
    
    with open(fname, 'wb') as output:
        pickle.dump(gobj, output, -1)

def load_gpcobj_pkl(fname):
    """ loading gpc object including infos about input pdfs, polynomials, grid etc...
    
    load_gpc(gobj, filename)
    
    fname ... filename with .pkl extension """
    
    with open(fname, 'rb') as input:
        return pickle.load(input)  

def save_data_txt(data, fname):
    """ saving data (quantity of interest) in .txt file (e.g. coeffs, mean, std, ...)
    
    save_data_txt(qoi, filename)
    
    data  ... quantity of interest, 2D np.array() 
    fname ... filename with .txt extension """
    
    np.savetxt(fname, data, fmt='%.10e', delimiter='\t', newline='\n', header='', footer='')

def load_data_hdf5(data, fname, loc):
    """ loading quantity of interest from .hdf5 file (e.g. coeffs, mean, std, ...)
    
    save_data_hdf5(data, filename, loc)
    
    data   ... quantity of interest, 2D np.array() 
    fname  ... filename with .hdf5 extension
    loc    ... location (folder and name) in hdf5 file e.g. data/phi """
    
    with h5py.File(fname, 'r') as f:
        d=f[loc]        
        return d
        
def save_data_hdf5(data, fname, loc):
    """ saving quantity of interest in .hdf5 file (e.g. coeffs, mean, std, ...)
    
    save_data_hdf5(data, filename, loc)
    
    data   ... quantity of interest, 2D np.array() 
    fname  ... filename with .hdf5 extension
    loc    ... location (folder and name) in hdf5 file e.g. data/phi """
    
    with h5py.File(fname, 'a') as f:
        f.create_dataset(loc, data=data)

def save_sobol_idx(sobol_idx, fname):
    """ saving sobol_idx list
    
    save_sobol_idx(sobol_idx, filename)
    
    sobol_idx ... sobol indices from method .sobol, list of nparray
    filename  ... filename with .txt extension """
    
    f = open(fname, 'w')
    f.write('# Parameter index list of Sobol indices:\n')
    for i_line in range(len(sobol_idx)):
        for i_entry in range(len(sobol_idx[i_line])):
            if i_entry > 0:
                f.write(', ')
            f.write('{}'.format(sobol_idx[i_line][i_entry]))
        if i_line < range(len(sobol_idx)):    
            f.write('\n')
        
    f.close()    
    
def read_sobol_idx(fname):
    """ reading sobol_idx list
    
    read_sobol_idx(filename)
    
    Input: 
        fname     ... filename with .txt extension
    
    Output:
        sobol_idx ... sobol indices, list of nparray
    """
     
    f = open(fname,'r')

    line = f.readline().strip('\n')
    sobol_idx = []

    while line:
        
        # ignore comments in textfile
        if line[0] == '#':
           line = f.readline().strip('\n')
           continue
        
        else:
            # read comma separated indices and convert to nparray
            sobol_idx.append(np.array([int(x) for x in line.split(',') if x]))
      
        line = f.readline().strip('\n')
    
    return sobol_idx    
    
def wrap_function(function, x, args):
    """ function wrapper to call anonymous function with variable number of arguments (tuple)
    
    save_qoi(qoi, filename)
    
    qoi   ... quantity of interest, 2D np.array() 
    fname ... filename with .txt extension"""
    def function_wrapper(*wrapper_args):
        return function(*(wrapper_args + x + args))
    
    return function_wrapper
            
def calc_Nc(order, dim):
    """ calc_Nc calculates the number of PCE coefficients by the used order
        and dimension.
      
        Nc   = (order+dim)! / (order! * dim!)
 
        Nc = calc_Nc( order , dim )
 
        Input:
            order   ... order of expansion
            dim     ... Number of random variables
 
        Output:
            Nc      ... Number of coefficients """
            
    return scipy.special.factorial(order+dim) / (scipy.special.factorial(order) * scipy.special.factorial(dim))

def calc_Nc_sparse(p_d, p_g, p_i, dim):
    """ calc_Nc_sparse calculates the number of PCE coefficients for a specific
        maximum order in each dimension p_d, maximum order of interacting 
        polynomials p_g and the interaction order p_i.
      
        Nc = calc_Nc_sparse(p_d, p_g, p_i, dim)
 
        Input:
            p_d ... maximum order in each dimension (scalar or np.array)
            p_g ... maximum order of interacting polynomials
            p_i ... interaction order            
            dim ... Number of dimensions
 
        Output:
            Nc      ... Number of coefficients """
    
    p_d = np.array(p_d)
    
    if p_d.size==1:
        p_d = p_d*np.ones(dim)
            
    # generate multi-index list up to maximum order
    if dim == 1:
        poly_idx = np.array([np.linspace(0,p_d,p_d+1)]).astype(int).transpose()
    else:
        poly_idx = allVL1leq(int(dim), p_g)
        
    for i_dim in range(dim):
        # add multi-indexes to list when not yet included
        if p_d[i_dim] > p_g:
            poly_add_dim = np.linspace(p_g+1, p_d[i_dim], p_d[i_dim]-(p_g+1) + 1)
            poly_add_all = np.zeros([poly_add_dim.shape[0],dim])
            poly_add_all[:,i_dim] = poly_add_dim               
            poly_idx = np.vstack([poly_idx,poly_add_all.astype(int)])
                
        # delete multi-indexes from list when they exceed individual max order of parameter     
        elif p_d[i_dim] < p_g:    
            poly_idx = poly_idx[poly_idx[:,i_dim]<=p_d[i_dim],:]
                
    # Consider interaction order (filter out multi-indices exceeding it)
    poly_idx = poly_idx[np.sum(poly_idx>0,axis=1)<=p_i,:]
        
    return poly_idx.shape[0]

def pdf_beta(x,p,q,a,b):
    """ pdf_beta calcualtes the probability density function of the beta 
        distribution in the inverval [a,b]

        pdf = pdf_beta( x , p , q , a , b )

        pdf = (gamma(p)*gamma(q)/gamma(p+q).*(b-a)**(p+q-1))**(-1) *
              (x-a)**(p-1) * (b-x)**(q-1);

        Input:
                x ... array of values of random variable x
                a ... MIN boundary 
                b ... MAX boundary
                p ... parameter defining the distribution shape
                q ... parameter defining the distribution shape

       Output:
               pdf ... value of probability density function """
    return (scipy.special.gamma(p)*scipy.special.gamma(q)/scipy.special.gamma(p+q)\
        * (b-a)**(p+q-1))**(-1) * (x-a)**(p-1) * (b-x)**(q-1)


def run_reg_adaptive(pdftype, pdfshape, limits, func, args=(), order_start=0, order_end=10, eps=1E-3, print_out=False):
    """  
    Adaptive regression approach based on leave one out cross validation error
    estimation
    
    Parameters
    ----------

    pdftype : list
              Type of probability density functions of input parameters,
              i.e. ["beta", "norm",...]
    pdfshape : list of lists
               Shape parameters of probability density functions
               s1=[...] "beta": p, "norm": mean
               s2=[...] "beta": q, "norm": std
               pdfshape = [s1,s2]
    limits : list of lists
             Upper and lower bounds of random variables (only "beta")
             a=[...] "beta": lower bound, "norm": n/a define 0
             b=[...] "beta": upper bound, "norm": n/a define 0
             limits = [a,b]
    func : callable func(x,*args)
           The objective function to be minimized.
    args : tuple, optional
           Extra arguments passed to func, i.e. f(x,*args).         
    order_start : int, optional
                  Initial gpc expansion order (maximum order)
    order_end : int, optional
                Maximum gpc expansion order to expand to
    eps : float, optional
          Relative mean error of leave one out cross validation
    print_out : boolean, optional
          Print output of iterations and subiterations (True/False)      

    Returns
    -------
    gobj : object
           gpc object
    res  : ndarray
           Funtion values at grid points of the N_out output variables
           size: [N_grid x N_out]        
    """
    
    # initialize iterators
    eps_gpc = eps+1
    i_grid = 0
    i_iter = 0
    interaction_order_count = 0
    interaction_order_max = -1
    DIM = len(pdftype)
    order = order_start
    
    while (eps_gpc > eps):
        
        # reset sub iteration if all polynomials of present order were added to gpc expansion
        if interaction_order_count > interaction_order_max:
            interaction_order_count = 0
            i_iter = i_iter + 1
            if print_out:    
                print("Iteration #{}".format(i_iter))
                print("=============\n")
           
        if i_iter == 1:
            # initialize gPC object in first iteration
            grid_init = randomgrid(pdftype, pdfshape, limits, np.ceil(1.2*calc_Nc(order,DIM)))
            regobj = reg(pdftype, pdfshape, limits, order*np.ones(DIM), order_max = order, interaction_order = DIM, grid = grid_init)
        else:
                       
            # generate new multi-indices if maximum gpc order increased
            if interaction_order_count == 0:
                order = order + 1
                if (order > order_end):
                    return regobj
                poly_idx_all_new = allVL1leq(DIM, order)
                poly_idx_all_new = poly_idx_all_new[np.sum(poly_idx_all_new,axis=1) == order]
                interaction_order_list = np.sum(poly_idx_all_new > 0,axis=1)
                interaction_order_max = np.max(interaction_order_list)
                interaction_order_count = 1
            
            if print_out: 
                print("   Subiteration #{}".format(interaction_order_count))
                print("   =============\n")
            # filter out polynomials of interaction_order = interaction_order_count
            poly_idx_added = poly_idx_all_new[interaction_order_list==interaction_order_count,:]
            
            # add polynomials to gpc expansion
            regobj.enrich_polynomial_basis(poly_idx_added)
            
            # generate new grid-points
            regobj.enrich_gpc_matrix_samples(1.2)
            
            interaction_order_count = interaction_order_count + 1    
            
        # run repeated simulations
        for i_grid in range(i_grid,regobj.grid.coords.shape[0]):
            if print_out:             
                print("   Performing simulation #{}\n".format(i_grid+1))
            # read conductivities from grid
            x = regobj.grid.coords[i_grid,:]
            
            # evaluate function at grid point
            res = func(x, *(args))
            
            # append result to solution matrix (RHS)
            if i_grid == 0:
                RES = res
            else:
                RES = np.vstack([RES, res])
        
        # increase grid counter by one for next iteration (to not repeat last simulation)         
        i_grid = i_grid + 1
        
        # perform leave one out cross validation 
        eps_gpc = regobj.LOOCV(RES)
        if print_out: 
            print("    -> relerror_LOOCV = {}\n".format(eps_gpc))
        
    return regobj, RES
    
#%%############################################################################
# gpc object class
###############################################################################

class gpc:
    def __init__(self):
        """ Initialize gpc class """
        dummy = 1
    
    def setup_polynomial_basis(self):
        """ Setup polynomial basis functions for a maximum order expansion """
        #print 'Setup polynomial basis functions ...'
        # Setup list of polynomials and their coefficients up to the desired order
        #
        #  poly    |     DIM_1     DIM_2    ...    DIM_M
        # -----------------------------------------------
        # Poly_1   |  [coeffs]  [coeffs]   ...  [coeffs]          
        # Poly_2   |  [coeffs]  [coeffs]   ...  [coeffs]
        #   ...    |  [coeffs]  [coeffs]   ...   [0]
        #   ...    |  [coeffs]  [coeffs]   ...   [0]
        #   ...    |  [coeffs]  [coeffs]   ...   ...
        # Poly_No  |   [0]      [coeffs]   ...   [0]
        #
        # size: [Max individual order x DIM]   (includes polynomials also not used)
        
        Nmax = int(np.max(self.order))
        
        # 2D list of polynomials (lookup)
        self.poly      = [[0 for x in range(self.DIM)] for x in range(Nmax+1)]
        
        # 2D array of polynomial normalization factors (lookup)
        # [Nmax+1 x DIM]
        self.poly_norm = np.zeros([Nmax+1,self.DIM])
        
        for i_DIM in range(self.DIM):
            for i_order in range(Nmax+1):
                if self.pdftype[i_DIM] == "beta": 
                    p = self.pdfshape[0][i_DIM] # beta-distr: alpha=p /// jacobi-poly: alpha=q-1  !!!
                    q = self.pdfshape[1][i_DIM] # beta-distr: beta=q  /// jacobi-poly: beta=p-1   !!!
                        
                    # determine polynomial normalization factor
                    beta_norm = (scipy.special.gamma(q)*scipy.special.gamma(p)/scipy.special.gamma(p+q)*(2.0)**(p+q-1))**(-1)
                    jacobi_norm = 2**(p+q-1) / (2.0*i_order+p+q-1)*scipy.special.gamma(i_order+p)*scipy.special.gamma(i_order+q) / (scipy.special.gamma(i_order+p+q-1)*scipy.special.factorial(i_order))
                    self.poly_norm[i_order,i_DIM] = (jacobi_norm * beta_norm)
                        
                    # add entry to polynomial lookup table
                    self.poly[i_order][i_DIM] = scipy.special.jacobi(i_order, q-1, p-1, monic=0)/np.sqrt(self.poly_norm[i_order,i_DIM]) 
                    
#==============================================================================
#                     alpha = self.pdfshape[0][i_DIM]-1 # here: alpha' = p-1
#                     beta = self.pdfshape[1][i_DIM]-1  # here: beta' = q-1
#                         
#                     # determine polynomial normalization factor
#                     beta_norm = (scipy.special.gamma(beta+1)*scipy.special.gamma(alpha+1)/scipy.special.gamma(beta+1+alpha+1)*(2.0)**(beta+1+alpha+1-1))**(-1)
#                     jacobi_norm = 2**(alpha+beta+1) / (2.0*i_order+alpha+beta+1)*scipy.special.gamma(i_order+alpha+1)*scipy.special.gamma(i_order+beta+1) / (scipy.special.gamma(i_order+alpha+beta+1)*scipy.special.factorial(i_order))
#                     self.poly_norm[i_order,i_DIM] = (jacobi_norm * beta_norm)
#                         
#                     # add entry to polynomial lookup table
#                     self.poly[i_order][i_DIM] = scipy.special.jacobi(i_order, beta, alpha, monic=0) # ! alpha beta changed definition in scipy!                  
#==============================================================================
                    
                if self.pdftype[i_DIM] == "normal" or self.pdftype[i_DIM] == "norm":
                        
                    # determine polynomial normalization factor
                    hermite_norm = scipy.special.factorial(i_order)
                    self.poly_norm[i_order,i_DIM] = hermite_norm
                        
                    # add entry to polynomial lookup table
                    self.poly[i_order][i_DIM] = scipy.special.hermitenorm(i_order, monic=0)/np.sqrt(self.poly_norm[i_order,i_DIM]) 
                        
        
        # Determine 2D multi-index array (order) of basis functions w.r.t. 2D array
        # of polynomials self.poly
        #
        # poly_idx |     DIM_1       DIM_2       ...    DIM_M
        # -------------------------------------------------------
        # basis_1  |  [order_D1]  [order_D2]     ...  [order_DM]    
        # basis_2  |  [order_D1]  [order_D2]     ...  [order_DM]
        #  ...     |  [order_D1]  [order_D2]     ...  [order_DM]
        #  ...     |  [order_D1]  [order_D2]     ...  [order_DM]
        #  ...     |  [order_D1]  [order_D2]     ...  [order_DM]
        # basis_Nb |  [order_D1]  [order_D2]     ...  [order_DM]
        #
        # size: [No. of basis functions x DIM]
        
        # generate multi-index list up to maximum order
        if self.DIM == 1:
            self.poly_idx = np.array([np.linspace(0,self.order_max,self.order_max+1)]).astype(int).transpose()
        else:
            self.poly_idx = allVL1leq(self.DIM, self.order_max)
            
        
        for i_DIM in range(self.DIM):
            # add multi-indexes to list when not yet included
            if self.order[i_DIM] > self.order_max:
                poly_add_dim = np.linspace(self.order_max+1, self.order[i_DIM], self.order[i_DIM]-(self.order_max+1) + 1)
                poly_add_all = np.zeros([poly_add_dim.shape[0],self.DIM])
                poly_add_all[:,i_DIM] = poly_add_dim               
                self.poly_idx = np.vstack([self.poly_idx,poly_add_all.astype(int)])
            # delete multi-indexes from list when they exceed individual max order of parameter     
            elif self.order[i_DIM] < self.order_max:    
                self.poly_idx = self.poly_idx[self.poly_idx[:,i_DIM]<=self.order[i_DIM],:]
                
        # Consider interaction order (filter out multi-indices exceeding it)
        if self.interaction_order < self.DIM: 
            self.poly_idx = self.poly_idx[np.sum(self.poly_idx>0,axis=1)<=self.interaction_order,:]        
        
        self.N_poly = self.poly_idx.shape[0]
           
#==============================================================================
#         x1 = [np.array(range(self.order[i]+1)) for i in range(self.DIM)]
#         self.poly_idx = []
#         
#         for element in itertools.product(*x1):
#             if np.sum(element) <= self.maxorder:
#                 self.poly_idx.append(element)
#                 
#         self.poly_idx = np.array(self.poly_idx) 
#==============================================================================
        
#==============================================================================
#         x1 = [np.array(range(self.order[i]+1)) for i in range(self.DIM)]
#         order_idx    = misc.combvec(x1)
#         
#         # filter for individual maximum expansion order
#         for i in range(self.DIM):
#             order_idx = order_idx[order_idx[:,i] <= self.order[i]]
# 
#         # filter for total maximum order
#         self.poly_idx = order_idx[np.sum(order_idx,axis = 1) <= self.order_max]
#==============================================================================
        
        # construct array of scaling factors to normalize basis functions <psi^2> = int(psi^2*p)dx
        # [Npolybasis x 1]
        self.poly_norm_basis = np.ones([self.poly_idx.shape[0],1])
        for i_poly in range(self.poly_idx.shape[0]):
            for i_DIM in range(self.DIM):
                self.poly_norm_basis[i_poly] *= self.poly_norm[self.poly_idx[i_poly,i_DIM],i_DIM]
    
    def enrich_polynomial_basis(self, poly_idx_added, form_A=True):
        """ Enrich polynomial basis functions and add new columns to gpc matrix 
            poly_idx_added ... array of added polynomials (order), np.array() [N_poly_added x DIM]"""
        
        # determine if polynomials in poly_idx_added are already present in self.poly_idx if so, delete them
        poly_idx_tmp = []
        for new_row in poly_idx_added:
            not_in_poly_idx = True
            for row in self.poly_idx:
                if np.allclose(row, new_row):
                    not_in_poly_idx = False
            if not_in_poly_idx:
                poly_idx_tmp.append(new_row)
        
        # if all polynomials are already present end routine
        if len(poly_idx_tmp) == 0:        
            return
        else:            
            poly_idx_added = np.vstack(poly_idx_tmp)
        
        # determine highest order added        
        order_max_added = np.max(np.max(poly_idx_added))
        
        # get current maximum order 
        order_max_current = len(self.poly)-1
        
        # Append list of polynomials and their coefficients up to the desired order
        #
        #  poly    |     DIM_1     DIM_2    ...    DIM_M
        # -----------------------------------------------
        # Poly_1   |  [coeffs_old]  [coeffs_old]   ...  [coeffs_old]          
        # Poly_2   |  [coeffs_old]  [coeffs_old]   ...  [coeffs_old]
        #   ...    |  [coeffs_old]  [coeffs_old]   ...  [coeffs_old]
        #   ...    |  [coeffs_old]  [coeffs_old]   ...  [coeffs_old]
        #   ...    |      ...           ...        ...      ...
        # Poly_No  |  [coeffs_new]  [coeffs_new]   ...  [coeffs_new]
        #
        # size: [Max order x DIM]
        
        # preallocate new rows to polynomial lists
        for i in range(order_max_added-order_max_current):
            self.poly.append([0 for x in range(self.DIM)])
            self.poly_norm = np.vstack([self.poly_norm, np.zeros(self.DIM)])
                
        for i_DIM in range(self.DIM):
            for i_order in range(order_max_current+1,order_max_added+1):
                if self.pdftype[i_DIM] == "beta":
                    p = self.pdfshape[0][i_DIM]
                    q = self.pdfshape[1][i_DIM]
                        
                    # determine polynomial normalization factor
                    beta_norm = (scipy.special.gamma(q)*scipy.special.gamma(p)/scipy.special.gamma(p+q)*(2.0)**(p+q-1))**(-1)
                    jacobi_norm = 2**(p+q-1) / (2.0*i_order+p+q-1)*scipy.special.gamma(i_order+p)*scipy.special.gamma(i_order+q) / (scipy.special.gamma(i_order+p+q-1)*scipy.special.factorial(i_order))
                    self.poly_norm[i_order,i_DIM] = (jacobi_norm * beta_norm)
                        
                    # add entry to polynomial lookup table
                    self.poly[i_order][i_DIM] = scipy.special.jacobi(i_order, q-1, p-1, monic=0)/np.sqrt(self.poly_norm[i_order,i_DIM])  # ! beta = p-1 and alpha=q-1 (consider definition in scipy.special.jacobi !!) 
                    
                if self.pdftype[i_DIM] == "normal" or self.pdftype[i_DIM] == "norm":
                        
                    # determine polynomial normalization factor
                    hermite_norm = scipy.special.factorial(i_order)
                    self.poly_norm[i_order,i_DIM] = hermite_norm
                        
                    # add entry to polynomial lookup table
                    self.poly[i_order][i_DIM] = scipy.special.hermitenorm(i_order, monic=0)/np.sqrt(self.poly_norm[i_order,i_DIM]) 
        
        # append new multi-indexes to old poly_idx array
        self.poly_idx = np.vstack([self.poly_idx, poly_idx_added])
        #self.poly_idx = unique_rows(self.poly_idx)
        self.N_poly = self.poly_idx.shape[0]
        
        # extend array of scaling factors to normalize basis functions <psi^2> = int(psi^2*p)dx
        # [Npolybasis x 1]
        N_poly_new = poly_idx_added.shape[0]
        poly_norm_basis_new = np.ones([N_poly_new,1])
        for i_poly in range(N_poly_new):
            for i_DIM in range(self.DIM):
                poly_norm_basis_new[i_poly] *= self.poly_norm[poly_idx_added[i_poly,i_DIM],i_DIM]
        
        self.poly_norm_basis = np.vstack([self.poly_norm_basis, poly_norm_basis_new])
        
        if form_A:
            # append new columns to gpc matrix [N_grid x N_poly_new]
            A_new_columns = np.zeros([self.N_grid,N_poly_new])
            for i_poly_new in range(N_poly_new):
                A1 = np.ones(self.N_grid)
                for i_DIM in range(self.DIM):
                    A1 *= self.poly[poly_idx_added[i_poly_new][i_DIM]][i_DIM](self.grid.coords_norm[:,i_DIM])
                A_new_columns[:,i_poly_new] = A1
            
            self.A = np.hstack([self.A, A_new_columns])
            self.Ainv = np.linalg.pinv(self.A)
    
    def enrich_gpc_matrix_samples(self, N_samples_N_poly_ratio):
        """ add sample points according to input pdfs to grid and enrich gpc matrix 
            such that the ratio of rows/columns is N_samples_N_poly_ratio """
        
        # Number of new grid points
        N_grid_new = int(np.ceil(N_samples_N_poly_ratio * self.A.shape[1] - self.A.shape[0]))
        
        if N_grid_new > 0:
            # Generate new grid points
            newgridpoints = randomgrid(self.pdftype, self.pdfshape, self.limits, N_grid_new)
            
            # append points to existing grid
            self.grid.coords = np.vstack([self.grid.coords, newgridpoints.coords])
            self.grid.coords_norm = np.vstack([self.grid.coords_norm, newgridpoints.coords_norm])
            self.N_grid = self.grid.coords.shape[0]
            
            # determine new row of gpc matrix
            a = np.zeros([N_grid_new, self.N_poly])
            for i_poly in range(self.N_poly):
                a1 = np.ones(N_grid_new)
                for i_DIM in range(self.DIM):
                    a1 *= self.poly[self.poly_idx[i_poly][i_DIM]][i_DIM](newgridpoints.coords_norm[:,i_DIM])
                a[:,i_poly] = a1
            
            # append new row to gpc matrix    
            self.A = np.vstack([self.A,a])
            
            # invert gpc matrix Ainv [N_basis x N_grid]
            self.Ainv  = np.linalg.pinv(self.A) 

        
    def construct_gpc_matrix(self):
        """ construct the gpc matrix A [N_grid x N_poly] """
        #print 'Constructing gPC matrix ...'
                
        self.A = np.zeros([self.N_grid,self.N_poly])
        
        for i_poly in range(self.N_poly):
            A1 = np.ones(self.N_grid)
            for i_DIM in range(self.DIM):
                A1 *= self.poly[self.poly_idx[i_poly][i_DIM]][i_DIM](self.grid.coords_norm[:,i_DIM])
            self.A[:,i_poly] = A1
        
        # invert gpc matrix Ainv [N_basis x N_grid]
        self.Ainv  = np.linalg.pinv(self.A)    
        


#%%############################################################################
# Postprocessing methods
###############################################################################
    def mean(self,coeffs):
        """ calculate the expected value 
            
            mean = calc_mean(self, coeffs)
            
            input:  coeffs ... gpc coefficients, np.array() [N_coeffs x N_out]            
            
            output: mean   ... mean, nparray [1 x N_out]
        """
            
        mean = coeffs[0,:]
        mean = mean[np.newaxis,:]
        return mean
        
    def std(self,coeffs):
        """ calculate the standard deviation  
            
            std = calc_std(self, coeffs)
            
            input:  coeffs ... gpc coefficients, np.array() [N_coeffs x N_out]            
            
            output: std    ... standard deviation, nparray [1 x N_out]
        """
            
        # return np.sqrt(np.sum(np.multiply(np.square(self.coeffs[1:,:]),self.poly_norm_basis[1:,:]),axis=0))
        std = np.sqrt(np.sum(np.square(coeffs[1:,:]),axis=0))
        std = std[np.newaxis,:]
        return std
        
    def MC_sampling(self, coeffs, N_samples, output_idx=[]):
        """ randomly sample the gpc expansion to determine output pdfs in specific points   
            
            samples_in, samples_out = MC_sampling(self, coeffs, N_samples, output_idx=[])
            
            input:  coeffs      ... gpc coefficients, np.array() [N_coeffs x N_out]
                    N_samples   ... number of random samples drawn from the respective input pdfs
                    output_idx  ... (optional) idx of output quantities to consider
                                    (if not, all outputs are considered), np.array() [1 x N_out] 
            
            output: samples_in  ... generated samples, nparray [N_samples x DIM]
                    samples_out ... gpc solutions, nparray [N_samples x N_out]
        """
        self.N_out = coeffs.shape[1]
             
        # if output index list is not provided, sample all gpc ouputs 
        if not output_idx:
            output_idx = np.linspace(0,self.N_out-1,self.N_out)
            output_idx = output_idx[np.newaxis,:]
            
        #np.random.seed()        
        
        # generate random samples for each random input variable [N_samples x DIM]
        samples_in = np.zeros([N_samples, self.DIM]) 
        for i_DIM in range(self.DIM):
            if self.pdftype[i_DIM] == "beta":
                samples_in[:,i_DIM] = (np.random.beta(self.pdfshape[0][i_DIM], self.pdfshape[1][i_DIM],[N_samples,1])*2.0 - 1)[:,0]
            if self.pdftype[i_DIM] == "norm" or self.pdftype[i_DIM] == "normal":
                samples_in[:,i_DIM] = (np.random.normal(0, 1, [N_samples,1]))[:,0]
        
        samples_out = self.evaluate(coeffs, samples_in, output_idx)
        return samples_in, samples_out
    
    def evaluate(self, coeffs, xi, output_idx=[]):
        """ calculate gpc approximation in points with output_idx and normalized 
        parameters xi (interval: [-1, 1]) 
        example: nparray = evaluate( [[xi_1_p1 ... xi_DIM_p1] , 
                                      [xi_1_p2 ... xi_DIM_p2]], nparray([[0,5,13]]) )
        
        calc_pdf = evaluate(self, coeffs, xi, output_idx)    
        
        input:  coeffs     ... gpc coefficients, np.array() [N_coeffs x N_out]
                xi         ... normalized coordinates of random variables, np.array() [N x DIM]
                output_idx ... (optional) idx of output quantities to consider
                               (if not provided, all outputs are considered), np.array() [1 x N_out] 
        
        output: y ... gpc approximation at normalized coordinates xi, np.array() [N_xi x N_out] 
        """
        
        self.N_out = coeffs.shape[1]
        self.N_poly = self.poly_idx.shape[0]
        
        # if point index list is not provided, evaluate over all points 
        if len(output_idx) == 0:
            output_idx = np.arange(self.N_out, dtype=int)
            output_idx = output_idx[np.newaxis,:]        
        
        N_out_eval = output_idx.shape[1]
        N_x = xi.shape[0]
        
        y = np.zeros([N_x, N_out_eval])
        for i_poly in range(self.N_poly):
            A1 = np.ones(N_x)
            for i_DIM in range(self.DIM):
                A1 *= self.poly[self.poly_idx[i_poly][i_DIM]][i_DIM](xi[:,i_DIM])
            y += np.outer(A1, coeffs[i_poly, output_idx.astype(int)])

        return y  
        
    def sobol(self, coeffs):
        """ Determine the available sobol indices 
        
        sobol, sobol_idx = calc_sobol(self, coeffs)    
        
        input:  coeffs    ... gpc coefficients, np.array() [N_coeffs x N_out]
                
        output: sobol     ... unnormalized sobol_indices, nparray [N_sobol x N_out] 
                sobol_idx ... sobol_idx list indicating the parameter combinations 
                              in rows of sobol, list of np.array [N_sobol x DIM]  
        """
                    
        N_sobol_theoretical = 2**self.DIM - 1
        N_coeffs = coeffs.shape[0]
        
        if N_coeffs == 1:
            raise Exception('Number of coefficients is 1 ... no sobol indices to calculate ...')
            
        # Generate boolean matrix of all basis functions where order > 0 = True
        # size: [N_coeffs x DIM] 
        sobol_mask = self.poly_idx != 0
        
        # look for unique combinations (i.e. available sobol combinations)
        # size: [N_sobol x DIM]
        sobol_idx_bool = unique_rows(sobol_mask)
        
        # delete the first row where all polys are order 0 (no sensitivity)
        sobol_idx_bool = np.delete(sobol_idx_bool,[0],axis=0)
        N_sobol_available = sobol_idx_bool.shape[0] 
        
        # check which basis functions contribute to which sobol coefficient set 
        # True for specific coeffs if it contributes to sobol coefficient
        # size: [N_coeffs x N_sobol]
        sobol_poly_idx = np.zeros([N_coeffs,N_sobol_available])
        for i_sobol in range(N_sobol_available):
            sobol_poly_idx[:,i_sobol] =  np.all(sobol_mask == sobol_idx_bool[i_sobol], axis=1)
            
        # calculate sobol coefficients matrix by summing up the individual
        # contributions to the respective sobol coefficients
        # size [N_sobol x N_points]    
        sobol = np.zeros([N_sobol_available,coeffs.shape[1]])        
        for i_sobol in range(N_sobol_available):
            sobol[i_sobol,:] = np.sum(np.square(coeffs[sobol_poly_idx[:,i_sobol]==1,:]),axis=0)
            # not normalized polynomials:             
            # sobol[i_sobol,:] = np.sum(np.multiply(np.square(coeffs[sobol_poly_idx[:,i_sobol]==1,:]),self.poly_norm_basis[sobol_poly_idx[:,i_sobol]==1,:]),axis=0)  
           
        # sort sobol coefficients in descending order (w.r.t. first output only ...)
        idx_sort_descend_1st = np.argsort(sobol[:,0],axis=0)[::-1]
        sobol = sobol[idx_sort_descend_1st,:]
        sobol_idx_bool = sobol_idx_bool[idx_sort_descend_1st]
        
        # get list of sobol indices
        sobol_idx = [0 for x in range(sobol_idx_bool.shape[0])]
        for i_sobol in range(sobol_idx_bool.shape[0]):      
            sobol_idx[i_sobol] = np.array([i for i, x in enumerate(sobol_idx_bool[i_sobol,:]) if x])
         
        
        return sobol, sobol_idx
        
    def globalsens(self, coeffs):
        """ Determine the global derivative based sensitivity coefficients
        
        Reference:
        D. Xiu, Fast Numerical Methods for Stochastic Computations: A Review, 
        Commun. Comput. Phys., 5 (2009), pp. 242-272 eq. (3.14) page 255
        
        globalsens = calc_globalsens(self, coeffs)    
        
        input:  coeffs     ... gpc coefficients, np.array() [N_coeffs x N_out]
        
        output: globalsens ... global sensitivity coefficients, nparray [DIM x N_out]    
        """
        
        Nmax = int(len(self.poly))
        
        self.poly_der = [[0 for x in range(self.DIM)] for x in range(Nmax)]
        poly_der_int = [[0 for x in range(self.DIM)] for x in range(Nmax)]
        poly_int = [[0 for x in range(self.DIM)] for x in range(Nmax)]
        knots_list_1D = [0 for x in range(self.DIM)]
        weights_list_1D = [0 for x in range(self.DIM)]
        
        # generate quadrature points for numerical integration for each random
        # variable separately
        
        for i_DIM in range(self.DIM):
            if self.pdftype[i_DIM] == 'beta':    # Jacobi polynomials
                knots_list_1D[i_DIM], weights_list_1D[i_DIM] = quadrature_jacobi_1D(Nmax,self.pdfshape[0][i_DIM]-1, self.pdfshape[1][i_DIM]-1)
            if self.pdftype[i_DIM] == 'norm' or self.pdftype[i_DIM] == "normal":   # Hermite polynomials
                knots_list_1D[i_DIM], weights_list_1D[i_DIM] = quadrature_hermite_1D(Nmax)
        
        # preprocess polynomials        
        for i_DIM in range(self.DIM):
            for i_order in range(Nmax):
                
                # evaluate the derivatives of the polynomials
                self.poly_der[i_order][i_DIM] = np.polyder(self.poly[i_order][i_DIM])
                
                # evaluate poly and poly_der at quadrature points and integrate w.r.t. pdf (multiply with weights and sum up)
                # saved like self.poly [N_order x DIM]
                poly_int[i_order][i_DIM]     = np.sum(np.dot(self.poly[i_order][i_DIM](knots_list_1D[i_DIM]), weights_list_1D[i_DIM]))
                poly_der_int[i_order][i_DIM] = np.sum(np.dot(self.poly_der[i_order][i_DIM](knots_list_1D[i_DIM]) , weights_list_1D[i_DIM]))
        
        N_poly = self.poly_idx.shape[0]
        poly_der_int_mult = np.zeros([self.DIM, N_poly])
        
        for i_sens in range(self.DIM):        
            for i_poly in range(N_poly):
                A1 = 1
                
                # evaluate complete integral                 
                for i_DIM in range(self.DIM):
                    if i_DIM == i_sens:
                        A1 *= poly_der_int[self.poly_idx[i_poly][i_DIM]][i_DIM]
                    else:
                        A1 *= poly_int[self.poly_idx[i_poly][i_DIM]][i_DIM]
                
                
                poly_der_int_mult[i_sens,i_poly] = A1
        
        # sum up over all coefficients        
        # [DIM x N_points]  = [DIM x N_poly] * [N_poly x N_points]
        globalsens = np.dot(poly_der_int_mult, coeffs)/(2**self.DIM)
        
        return globalsens
        
    def localsens(self, coeffs, xi):
        """ Determine the local derivative based sensitivity coefficients in the
        point of operation xi (normalized coordinates!).
        
        example: xi = np.array([[0,0,...,0]]) size: [1 x DIM]
        
        localsens = calc_localsens(self, coeffs, xi)    
        
        input:  coeffs     ... gpc coefficients, np.array() [N_coeffs x N_out]
                xi         ... point in variable space to evaluate local sensitivity in 
                               (norm. coordinates) np.array() [1 x DIM]
        
        output: localsens ... local sensitivity coefficients, np.array() [DIM x N_out]  
        """
        Nmax = len(self.poly)
        
        self.poly_der = [[0 for x in range(self.DIM)] for x in range(Nmax+1)]
        poly_der_xi = [[0 for x in range(self.DIM)] for x in range(Nmax+1)]
        poly_opvals = [[0 for x in range(self.DIM)] for x in range(Nmax+1)]
        
        # preprocess polynomials        
        for i_DIM in range(self.DIM):
            for i_order in range(Nmax+1):
                
                # evaluate the derivatives of the polynomials
                self.poly_der[i_order][i_DIM] = np.polyder(self.poly[i_order][i_DIM])
                
                # evaluate poly and poly_der at point of operation
                poly_opvals[i_order][i_DIM] =  self.poly[i_order][i_DIM](xi[1,i_DIM])
                poly_der_xi[i_order][i_DIM] =  self.poly_der[i_order][i_DIM](xi[1,i_DIM])
        
        N_vals = 1
        poly_sens = np.zeros([self.DIM, self.N_poly])
        
        for i_sens in range(self.DIM):        
            for i_poly in range(self.N_poly):
                A1 = np.ones(N_vals)
                
                # construct polynomial basis according to partial derivatives                
                for i_DIM in range(self.DIM):
                    if i_DIM == i_sens:
                        A1 *= poly_der_xi[self.poly_idx[i_poly][i_DIM]][i_DIM]
                    else:
                        A1 *= poly_opvals[self.poly_idx[i_poly][i_DIM]][i_DIM]
                poly_sens[i_sens,i_poly] = A1
        
        # sum up over all coefficients        
        # [DIM x N_points]  = [DIM x N_poly]  *   [N_poly x N_points]
        localsens = np.dot(poly_sens,coeffs)   
        
        return localsens
    
    def pdf(self, coeffs, N_samples, output_idx=[]):
        """ Determine the estimated pdfs of the output quantities 
        
        calc_pdf = calc_pdf(self, coeffs, N_samples, output_idx)    
        
        input:  coeffs     ... gpc coefficients, np.array() [N_coeffs x N_out]
                N_samples  ... Number of samples used to estimate output pdf 
                output_idx ... (optional) idx of output quantities to consider
                               (if not, all outputs are considered), np.array() [1 x N_out] 
        
        output: pdf_x ... x-coordinates of output pdf (output quantity), nparray [100 x N_out]
                pdf_y ... y-coordinates of output pdf (probability density of output quantity), nparray [100 x N_out]  
        """
        
        self.N_out = coeffs.shape[1]
        
        # if output index array is not provided, determine pdfs of all outputs 
        if not output_idx:
            output_idx = np.linspace(0,self.N_out-1,self.N_out)
            output_idx = output_idx[np.newaxis,:]
    
        # sample gPC expansion
        samples_in, samples_out = self.MC_sampling(coeffs, N_samples, output_idx)
        
        # determine kernel density estimates using Gaussian kernel
        pdf_x = np.zeros([100,self.N_out])
        pdf_y = np.zeros([100,self.N_out])
        
        for i_out in range(coeffs.shape[1]):
            kde = scipy.stats.gaussian_kde(samples_out.transpose(), bw_method=0.1/samples_out[:,i_out].std(ddof=1))
            pdf_x[:,i_out] = np.linspace(samples_out[:,i_out].min(), samples_out[:,i_out].max(), 100)
            pdf_y[:,i_out] = kde(pdf_x[:,i_out])
            
        return pdf_x, pdf_y
        

        
#%%############################################################################
# Regression based gpc object subclass
###############################################################################

class reg(gpc):
    def __init__(self, pdftype, pdfshape, limits, order, order_max, interaction_order, grid):
        """
        regression gpc class
        -----------------------
        reg(pdftype, pdfshape, limits, order, order_max, interaction_order, grid)
        
        input:
            pdftype   ... type of pdf 'beta' or 'jacobi' list [1 x DIM]
            pdfshape  ... shape parameters of pdfs list [2 x DIM]
                            beta-dist:   [[alpha], [beta]    ]
                            normal-dist: [[mean],  [variance]] 

            limits    ... upper and lower bounds of random variables list [2 x DIM]
                            beta-dist:   [[a1 ...], [b1 ...]]
                            normal-dist: [[0 ... ], [0 ... ]] (not used)
            order     ... maximum individual expansion order list [1 x DIM]
                          generates individual polynomials also if maximum expansion order 
                          in order_max is exceeded    
            order_max ... maximum expansion order (sum of all exponents) [1]
                          the maximum expansion order considers the sum of the 
                          orders of combined polynomials only 
            interaction_order ... number of random variables, which can interact with each other [1]
                          all polynomials are ignored, which have an interaction order greater than the specified 
            grid      ... grid object gnerated in .grid.py including grid.coords and grid.coords_norm                                      
        """
        gpc.__init__(self)
        self.pdftype        = pdftype
        self.pdfshape       = pdfshape
        self.limits         = limits 
        self.order          = order                     
        self.order_max      = order_max
        self.interaction_order = interaction_order
        self.DIM            = len(pdftype)
        self.grid           = grid        
        self.N_grid         = grid.coords.shape[0]
        
        # setup polynomial basis functions
        self.setup_polynomial_basis()    
        
        # construct gpc matrix [Ngrid x Npolybasis]
        self.construct_gpc_matrix()
 
    def expand(self, data):
        """ Determine the gPC coefficients by the regression method
        input:    data ... results from simulations with N_out output quantities, np.array() [N_grid x N_out]
        output: coeffs ... gPC coefficients, np.array() [N_coeffs x N_out]
        """
        #print 'Determine gPC coefficients ...'
        self.N_out = data.shape[1]
        
        if data.shape[0] != self.Ainv.shape[1]:
            if data.shape[1] != self.Ainv.shape[1]:
                print("Please check format of input data: matrix [N_grid x N_out] !")
            else:
                data = data.T
        # coeffs    ... [N_coeffs x N_points] 
        # Ainv      ... [N_coeffs x N_grid]
        # data      ... [N_grid   x N_points]        
        return np.dot(self.Ainv,data)


    def LOOCV(self, data):
         """ Perform leave one out cross validation 
         input:    data ... results from simulations with N_out output quantities, np.array() [N_grid x N_out]
         output: relerror_LOOCV ... relative mean error of leave one out cross validation
         """
         
         self.relerror = np.zeros(data.shape[0])
         
         for i in range(data.shape[0]):
             # get mask of eliminated row                
             mask = np.arange(data.shape[0]) != i
             
             # invert reduced gpc matrix
             Ainv_LOO = np.linalg.pinv(self.A[mask,:])
             
             # determine gpc coefficients
             coeffs_LOO = np.dot(Ainv_LOO, data[mask,:])
             
             self.relerror[i] = np.linalg.norm(data[i,:] - np.dot(self.A[i,:], coeffs_LOO)) / np.linalg.norm(data[i,:])
             
         self.relerror_LOOCV = np.mean(self.relerror)
         
         return self.relerror_LOOCV
             
#%%############################################################################
# Quadrature based gpc object subclass
###############################################################################

class quad(gpc):
    def __init__(self, pdftype, pdfshape, limits, order, order_max, interaction_order, grid):
        gpc.__init__(self)
        self.pdftype        = pdftype
        self.pdfshape       = pdfshape
        self.limits         = limits
        self.order          = order
        self.order_max      = order_max
        self.interaction_order = interaction_order
        self.DIM            = len(pdftype)
        self.grid           = grid  
        self.N_grid         = grid.coords.shape[0]
        
        # setup polynomial basis functions
        self.setup_polynomial_basis()
        
        # construct gpc matrix [Ngrid x Npolybasis]
        self.construct_gpc_matrix()
           
    def expand(self, data):
        """ Determine the gPC coefficients by numerical integration
        input:    data ... results from simulations with N_out output quantities, np.array() [N_grid x N_out]
        output: coeffs ... gPC coefficients, np.array() [N_coeffs x N_out]
        """
        #print 'Determine gPC coefficients ...'
        self.N_out = data.shape[1]
        
        # check if quadrature rule (grid) fits to the distribution (pdfs)
        grid_pdf_fit = 1
        for i_DIM in range(self.DIM):
            if self.pdftype[i_DIM] == 'beta':
                if not (self.grid.gridtype[i_DIM] == 'jacobi'):
                    grid_pdf_fit = 0
                    break
            elif (self.pdftype[i_DIM] == 'norm') or (self.pdftype[i_DIM] == 'normal'):
                if not (self.grid.gridtype[i_DIM] == 'hermite'):
                    grid_pdf_fit = 0
                    break
    
        # if not, calculate joint pdf
        if not(grid_pdf_fit):
            joint_pdf = np.ones(self.grid.coords_norm.shape)
            
            for i_DIM in range(self.DIM):
                if self.pdftype[i_DIM] == 'beta':
                    joint_pdf[:,i_DIM] = pdf_beta(self.grid.coords_norm[:,i_DIM],self.pdfshape[0][i_DIM],self.pdfshape[1][i_DIM],-1,1)
                if (self.pdftype[i_DIM] == 'norm' or self.pdftype[i_DIM] == 'normal'):
                    joint_pdf[:,i_DIM] = scipy.stats.norm.pdf(self.grid.coords_norm[:,i_DIM])
            
            joint_pdf = np.array([np.prod(joint_pdf,axis=1)]).transpose()
            
            # weight data with the joint pdf
            data = data*joint_pdf*2**self.DIM
        
        # scale rows of gpc matrix with quadrature weights
        A_weighted = np.dot(np.diag(self.grid.weights), self.A)
        # scale = np.outer(self.grid.weights, 1./self.poly_norm_basis)        
        # A_weighted = np.multiply(self.A, scale)
        
        # determine gpc coefficients [N_coeffs x N_output]
        return np.dot(data.transpose(), A_weighted).transpose()           
        
