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
from scipy.fftpack import ifft
from .misc import combvec
from .misc import allVL1leq

def quadrature_jacobi_1D(N, b, a):
    """ Get knots and weights of Jacobi polynomials (beta distribution) """
    # make array to count N: 0, 1, ..., N-1
    N_arr = np.array(range(1, np.int(N), 1))
                
    # compose diagonals for companion matrix
    t01 = 1.0*(b-a) / (2+a+b)             
    t02 = 1.0*((b-a)*(a+b)) / ((2*N_arr+a+b)*(2*N_arr+2+a+b))
    t1  = np.append(t01, t02)
    t2  = np.sqrt((4.0*N_arr*(N_arr+a)*(N_arr+b)*(N_arr+a+b))/((2*N_arr-1+a+b)*(2*N_arr+a+b)**2*(2*N_arr+1+a+b)))    
            
    # compose companion matrix
    T = np.diag(t1) + np.diag(t2,1) + np.diag(t2,-1)

    # evaluate roots of polynomials (the abscissas are the roots of the
    # characteristic polynomial, i.d. the eigenvalues of the companion matrix)
    # the weights can be derived from the corresponding eigenvectors.
    eigvals, eigvecs = np.linalg.eig(T)
    idx_sorted       = np.argsort(eigvals)
    eigvals_sorted   = eigvals[idx_sorted]
        
    weights = 2.0*eigvecs[0,idx_sorted]**2
    knots   = eigvals_sorted
          
    return knots, weights
        
def quadrature_hermite_1D(N):
    """ Get knots and weights of prob. Hermite polynomials (normal distribution) """
    N = np.int(N)
    knots, weights = np.polynomial.hermite_e.hermegauss(N)
    weights = np.array(list(2.0*weights/np.sum(weights)))        
            
    return knots, weights
   
def quadrature_cc_1D(N):
    """ Computes the Clenshaw Curtis nodes and weights """
    N = np.int(N)        
    if N == 1:
        knots = 0
        weights = 2
    else:
        n = N - 1
        C = np.zeros((N,2))
        k = 2*(1+np.arange(np.floor(n/2)))
        C[::2,0] = 2/np.hstack((1, 1-k*k))
        C[1,1] = -n
        V = np.vstack((C,np.flipud(C[1:n,:])))
        F = np.real(ifft(V, n=None, axis=0))
        knots = F[0:N,1]
        weights = np.hstack((F[0,0],2*F[1:n,0],F[n,0]))
            
    return knots, weights

def quadrature_fejer1_1D(N):
    """ Computes the Fejer type 1 nodes and weights
    
        This method uses a direct approach.  The paper by Waldvogel
        exhibits a more efficient approach using Fourier transforms.
     
        Reference:
        Philip Davis, Philip Rabinowitz,
        Methods of Numerical Integration,
        Second Edition,
        Dover, 2007,
        ISBN: 0486453391 Titel anhand dieser ISBN in Citavi-Projekt übernehmen,
        LC: QA299.3.D28.

        Walter Gautschi,
        Numerical Quadrature in the Presence of a Singularity,
        SIAM Journal on Numerical Analysis,
        Volume 4, Number 3, 1967, pages 357-362.

        Joerg Waldvogel,
        Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
        BIT Numerical Mathematics,
        Volume 43, Number 1, 2003, pages 1-18.
   """
    N = np.int(N)
        
    theta = np.zeros ( N )

    for i in range ( 0, N ):
        theta[i] = float ( 2 * N - 1 - 2 * i ) * np.pi / float ( 2 * N )

    knots = np.zeros ( N )
        
    for i in range ( 0, N ):
        knots[i] = np.cos ( theta[i] )

    weights = np.zeros ( N )

    for i in range ( 0, N ):
        weights[i] = 1.0
        jhi = ( N // 2 )
        for j in range ( 0, jhi ):
            angle = 2.0 * float ( j + 1 ) * theta[i]
            weights[i] = weights[i] - 2.0 * np.cos ( angle ) / float ( 4 * ( j + 1 ) ** 2 - 1 )

    for i in range ( 0, N ):
        weights[i] = 2.0 * weights[i] / float ( N )

    return knots, weights

def quadrature_fejer2_1D(N):
    """ Computes the Fejer type 2 nodes and weights (Clenshaw Curtis without boundary nodes)
        
        This method uses a direct approach.  The paper by Waldvogel
        exhibits a more efficient approach using Fourier transforms.
        
        Reference:
        Philip Davis, Philip Rabinowitz,
        Methods of Numerical Integration,
        Second Edition,
        Dover, 2007,
        ISBN: 0486453391 Titel anhand dieser ISBN in Citavi-Projekt übernehmen,
        LC: QA299.3.D28.

        Walter Gautschi,
        Numerical Quadrature in the Presence of a Singularity,
        SIAM Journal on Numerical Analysis,
        Volume 4, Number 3, 1967, pages 357-362.

        Joerg Waldvogel,
        Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
        BIT Numerical Mathematics,
        Volume 43, Number 1, 2003, pages 1-18.
    """
    N = np.int(N)
        
    if ( N == 1 ):

        knots = np.array ( [ 0.0 ] )
        weights = np.array ( [ 2.0 ] )

    elif ( N == 2 ):

        knots = np.array ( [ -0.5, +0.5 ] )
        weights = np.array ( [  1.0,  1.0 ] )

    else:
            
        theta = np.zeros ( N )

        for i in range ( 0, N ):
            theta[i] = float ( N - i ) * np.pi / float ( N + 1 )

        knots = np.zeros ( N )
            
        for i in range ( 0, N ):
            knots[i] = np.cos ( theta[i] )

        weights = np.zeros ( N )

        for i in range ( 0, N ):

            weights[i] = 1.0
            jhi = ( ( N - 1 ) // 2 )
                
            for j in range ( 0, jhi ):
                angle = 2.0 * float ( j + 1 ) * theta[i]
                weights[i] = weights[i] - 2.0 * np.cos ( angle ) / float ( 4 * ( j + 1 ) ** 2 - 1 )
                p = 2 *  ( ( N + 1 ) // 2 ) - 1

            weights[i] = weights[i] - np.cos ( float ( p + 1 ) * theta[i] ) / float ( p )

        for i in range ( 0, N ):
            weights[i] = 2.0 * weights[i] / float ( N + 1 )

    return knots, weights
    

def quadrature_patterson_1D(N):
	""" Computes the nested Gauss-Patterson nodes and weights for N = 1,3,7,15,31 """
	x = np.zeros(N)
	w = np.zeros(N)
	
	if N == 1:

		x = 0.0;

		w = 2.0;

	elif N == 3:
		
		x[0] = -0.77459666924148337704;
		x[1] =  0.0;
		x[2] =  0.77459666924148337704;
		
		w[0] = 0.555555555555555555556;
		w[1] = 0.888888888888888888889;
		w[2] = 0.555555555555555555556;

	elif N == 7:

		x[0] = -0.96049126870802028342;
		x[1] = -0.77459666924148337704;
		x[2] = -0.43424374934680255800;
		x[3] =  0.0;
		x[4] =  0.43424374934680255800;
		x[5] =  0.77459666924148337704;
		x[6] =  0.96049126870802028342;

		w[0] = 0.104656226026467265194;
		w[1] = 0.268488089868333440729;
		w[2] = 0.401397414775962222905;
		w[3] = 0.450916538658474142345;
		w[4] = 0.401397414775962222905;
		w[5] = 0.268488089868333440729;
		w[6] = 0.104656226026467265194;

	elif N == 15:

		x[ 0] = -0.99383196321275502221;
		x[ 1] = -0.96049126870802028342;
		x[ 2] = -0.88845923287225699889;
		x[ 3] = -0.77459666924148337704;
		x[ 4] = -0.62110294673722640294;
		x[ 5] = -0.43424374934680255800;
		x[ 6] = -0.22338668642896688163;
		x[ 7] =  0.0;
		x[ 8] =  0.22338668642896688163;
		x[ 9] =  0.43424374934680255800;
		x[10] =  0.62110294673722640294;
		x[11] =  0.77459666924148337704;
		x[12] =  0.88845923287225699889;
		x[13] =  0.96049126870802028342;
		x[14] =  0.99383196321275502221;

		w[ 0] = 0.0170017196299402603390;
		w[ 1] = 0.0516032829970797396969;
		w[ 2] = 0.0929271953151245376859;
		w[ 3] = 0.134415255243784220360;
		w[ 4] = 0.171511909136391380787;
		w[ 5] = 0.200628529376989021034;
		w[ 6] = 0.219156858401587496404;
		w[ 7] = 0.225510499798206687386;
		w[ 8] = 0.219156858401587496404;
		w[ 9] = 0.200628529376989021034;
		w[10] = 0.171511909136391380787;
		w[11] = 0.134415255243784220360;
		w[12] = 0.0929271953151245376859;
		w[13] = 0.0516032829970797396969;
		w[14] = 0.0170017196299402603390;
		
	elif N == 31:
		
		x[ 0] = -0.99909812496766759766;
		x[ 1] = -0.99383196321275502221;
		x[ 2] = -0.98153114955374010687;
		x[ 3] = -0.96049126870802028342;
		x[ 4] = -0.92965485742974005667;
		x[ 5] = -0.88845923287225699889;
		x[ 6] = -0.83672593816886873550;
		x[ 7] = -0.77459666924148337704;
		x[ 8] = -0.70249620649152707861;
		x[ 9] = -0.62110294673722640294;
		x[10] = -0.53131974364437562397;
		x[11] = -0.43424374934680255800;
		x[12] = -0.33113539325797683309;
		x[13] = -0.22338668642896688163;
		x[14] = -0.11248894313318662575;
		x[15] =  0.0;
		x[16] =  0.11248894313318662575;
		x[17] =  0.22338668642896688163;
		x[18] =  0.33113539325797683309;
		x[19] =  0.43424374934680255800;
		x[20] =  0.53131974364437562397;
		x[21] =  0.62110294673722640294;
		x[22] =  0.70249620649152707861;
		x[23] =  0.77459666924148337704;
		x[24] =  0.83672593816886873550;
		x[25] =  0.88845923287225699889;
		x[26] =  0.92965485742974005667;
		x[27] =  0.96049126870802028342;
		x[28] =  0.98153114955374010687;
		x[29] =  0.99383196321275502221;
		x[30] =  0.99909812496766759766;
		
		w[ 0] = 0.00254478079156187441540;
		w[ 1] = 0.00843456573932110624631;
		w[ 2] = 0.0164460498543878109338;
		w[ 3] = 0.0258075980961766535646;
		w[ 4] = 0.0359571033071293220968;
		w[ 5] = 0.0464628932617579865414;
		w[ 6] = 0.0569795094941233574122;
		w[ 7] = 0.0672077542959907035404;
		w[ 8] = 0.0768796204990035310427;
		w[ 9] = 0.0857559200499903511542;
		w[10] = 0.0936271099812644736167;
		w[11] = 0.100314278611795578771;
		w[12] = 0.105669893580234809744;
		w[13] = 0.109578421055924638237;
		w[14] = 0.111956873020953456880;
		w[15] = 0.112755256720768691607;
		w[16] = 0.111956873020953456880;
		w[17] = 0.109578421055924638237;
		w[18] = 0.105669893580234809744;
		w[19] = 0.100314278611795578771;
		w[20] = 0.0936271099812644736167;
		w[21] = 0.0857559200499903511542;
		w[22] = 0.0768796204990035310427;
		w[23] = 0.0672077542959907035404;
		w[24] = 0.0569795094941233574122;
		w[25] = 0.0464628932617579865414;
		w[26] = 0.0359571033071293220968;
		w[27] = 0.0258075980961766535646;
		w[28] = 0.0164460498543878109338;
		w[29] = 0.00843456573932110624631;
		w[30] = 0.00254478079156187441540;
	else:
		print("Number of points does not match Gauss-Patterson quadrature rule ...")
		
	return x,w

def denorm(coords_norm,pdftype,gridshape,limits):
    """ Denormalize grid from standardized ([-1, 1] except hermite) to original parameter space for simulations """
    coords = np.zeros(coords_norm.shape)
    
    for i_DIM in range(coords_norm.shape[1]):
        #if gridtype[i_DIM] == 'jacobi' or gridtype[i_DIM] == 'cc' or gridtype[i_DIM] == 'fejer2':
        if pdftype[i_DIM] == "beta":
            coords[:,i_DIM] = (coords_norm[:,i_DIM] + 1) / 2 * (limits[1][i_DIM]-limits[0][i_DIM]) + limits[0][i_DIM]    
        #if gridtype[i_DIM] == 'hermite':
        if pdftype[i_DIM] == "norm" or pdftype[i_DIM] == "normal":   
            coords[:,i_DIM] = coords_norm[:,i_DIM]*gridshape[1][i_DIM] + gridshape[0][i_DIM]
    
    return coords
        

def norm(coords, pdftype, gridshape, limits):
    """ Normalize grid from original parameter (except hermite) to standardized ([-1, 1] space for simulations """
    coords_norm = np.zeros(coords.shape)

    for i_DIM in range(coords.shape[1]):
        # if gridtype[i_DIM] == 'jacobi' or gridtype[i_DIM] == 'cc' or gridtype[i_DIM] == 'fejer2':
        if pdftype[i_DIM] == "beta":
            coords_norm[:, i_DIM] = (coords[:, i_DIM] - limits[0][i_DIM])
            coords_norm[:, i_DIM] = coords_norm[:, i_DIM] / (limits[1][i_DIM]-limits[0][i_DIM])*2.0-1

            # if gridtype[i_DIM] == 'hermite':
        if pdftype[i_DIM] == "norm" or pdftype[i_DIM] == "normal":
            coords_norm[:, i_DIM] = (coords[:, i_DIM] - gridshape[0][i_DIM])/ gridshape[1][i_DIM] 
    return coords_norm

#%% full-grid object subclass   
###############################################################################
    
class tensgrid():
    def __init__(self, pdftype, gridtype, gridshape, limits, N):
        #print '-> Generate full tensored grid'  
        self.pdftype   = pdftype         # 'beta', 'normal'
        self.gridtype  = gridtype        # 'jacobi', 'hermite', 'cc', 'fejer2'
        self.gridshape = gridshape       # pdfshape: jacobi: -> [alpha and beta] hermite: -> [mean, variance]
        self.limits    = limits          # limits: [min, max]        
        self.N         = N               # number of nodes in each dimension
        self.DIM       = len(self.N)     # number of dimension
                
        # get knots and weights of polynomials in each dimension
        self.knots_DIM = []
        self.weights_DIM = []
        for i_DIM in range(self.DIM):
            if self.gridtype[i_DIM] == 'jacobi':    # jacobi polynomials
                knots_temp, weights_temp = quadrature_jacobi_1D(self.N[i_DIM], self.gridshape[0][i_DIM]-1, self.gridshape[1][i_DIM]-1)                
            if self.gridtype[i_DIM] == 'hermite':   # hermite polynomials
                knots_temp, weights_temp = quadrature_hermite_1D(self.N[i_DIM])
            if self.gridtype[i_DIM] == 'cc':        # Clenshaw Curtis
                knots_temp, weights_temp = quadrature_cc_1D(self.N[i_DIM])
            if self.gridtype[i_DIM] == 'fejer2':    # Fejer type 2 (Clenshaw Curtis without boundary nodes)
                knots_temp, weights_temp = quadrature_fejer2_1D(self.N[i_DIM])
            if self.gridtype[i_DIM] == 'patterson': # Gauss-Patterson (Nested Legendre rule)
                knots_temp, weights_temp = quadrature_patterson_1D(self.N[i_DIM])    
            
            self.knots_DIM.append(knots_temp)
            self.weights_DIM.append(weights_temp)
                    
        # combine coordinates to full tensored grid (all combinations)
        self.coords_norm = combvec(self.knots_DIM)
        
        # rescale normalized coordinates in case of normal distributions and "fejer2" or "cc" grids
        # +- 0.675 * sigma -> 50%
        # +- 1.645 * sigma -> 90%
        # +- 1.960 * sigma -> 95%        
        # +- 2.576 * sigma -> 99%
        # +- 3.000 * sigma -> 99.73%
        for i_DIM in range(self.DIM):
            if (self.pdftype[i_DIM] == "norm" or self.pdftype[i_DIM] == "normal") and (not(self.gridtype[i_DIM] == "hermite")): 
                self.coords_norm[:,i_DIM] = self.coords_norm[:,i_DIM]*1.960
        
        # determine combined weights of Gauss quadrature
        self.weights = np.prod(combvec(self.weights_DIM),axis=1)/(2.0**self.DIM)
        
        # denormalize grid to original parameter space
        self.coords = denorm(self.coords_norm,self.pdftype,self.gridshape,self.limits)
        
#%% sparse-grid object subclass   
###############################################################################
    
class sparsegrid():
    def __init__(self, pdftype, gridtype, gridshape, limits, level, level_max, interaction_order, order_sequence_type):
        #print '-> Generate sparse grid' 
        self.pdftype   = pdftype         # 'beta', 'normal'
        self.gridtype  = gridtype        # 'jacobi', 'hermite', 'cc', 'fejer2'
        self.gridshape = gridshape       # pdfshape: jacobi: -> [alpha and beta] hermite: -> [mean, variance]
        self.limits    = limits          # limits: [min, max]        
        self.level     = level           # number of levels in each dimension
        self.level_max = level_max       # global combined level maximum
        self.interaction_order = interaction_order
        self.order_sequence_type = order_sequence_type # 'lin', 'exp' type of order sequence (common: 'exp')        
        self.DIM       = len(self.level) # number of dimension
        
        # make multi-index list
        print("    Generating multi-indices ...")
        order_sequence = []
        level_sequence = []
        for i_dim in range(self.DIM):
            if self.gridtype[i_dim] == 'fejer2':
                level_sequence.append(range(1,self.level[i_dim]+1))
            else:
                level_sequence.append(range(self.level[i_dim]+1))
            if self.order_sequence_type == 'exp':       # order = 2**level + 1
                if self.gridtype[i_dim] == 'fejer2':    # start with order = 1 @ level = 1
                    order_sequence.append(2**(np.linspace(1,self.level[i_dim],self.level[i_dim]))-1)
                    order_sequence[i_dim][0] = 1
                elif self.gridtype[i_dim] == 'patterson': # start with order = 1 @ level = 0 [1,3,7,15,31,...]
                    order_sequence.append(2**(np.linspace(0,self.level[i_dim],self.level[i_dim]+1)+1)-1)
                else:                                   # start with order = 1 @ level = 0
                    order_sequence.append(2**np.linspace(0,self.level[i_dim],self.level[i_dim]+1)+1)
                    order_sequence[i_dim][0] = 1
            elif self.order_sequence_type == 'lin':     # order = level
                if self.gridtype[i_dim] == 'fejer2':    # start with level = 1 @ order = 1           
                    order_sequence.append(np.linspace(1,self.level[i_dim]+1,self.level[i_dim]+1))
                elif self.gridtype[i_dim] == 'patterson': # start with order = 1 @ level = 0 [1,3,7,15,31,...]
                    print("Not possible in case of Gauss-Patterson grid ...")
                else:                                   # start with 
                    order_sequence.append(np.linspace(1,self.level[i_dim]+1,self.level[i_dim]+1))

#==============================================================================
#         l_order = misc.combvec(order_sequence)
#         l_level = misc.combvec(level_sequence)
# 
#         # filter multi-index list to max. dimension of according level and global maxlevel
#         # count number of Fejer dimensions
#         N_fejer = 0
#         for i_dim in range(self.DIM):
#             if self.gridtype[i_dim] == 'fejer2':
#                 N_fejer += 1
#         
#         maxlevel_dim_idx = []
#         for i_dim in range(self.DIM):
#             if not self.gridtype[i_dim]=='fejer2':
#                 offset = N_fejer
#             else:
#                 offset = N_fejer-1
#             maxlevel_dim_idx.append(np.array(np.sum(l_level,axis=1)==(l_level[:,i_dim]+offset)))   
# 
#         # delete entries exceeding global level maximum but keep individual max levels
#         maxlevel_global_idx = np.sum(l_level,axis=1)<=(self.level_max )
# 
#         filter_idx = maxlevel_global_idx
# 
#         for i_dim in range(self.DIM):
#             filter_idx |= maxlevel_dim_idx[i_dim] 
# 
#         l_order = l_order[filter_idx,:]
#         l_level = l_level[filter_idx,:]
#         
#         self.l_level = l_level
#         self.l_order = l_order
#==============================================================================
        
        if any("fejer2" in s for s in self.gridtype):
            if self.DIM == 1:
                l_level = np.array([np.linspace(1,self.level_max,self.level_max)]).transpose()
            else:    
                l_level = allVL1leq(self.DIM, self.level_max-self.DIM)
                l_level = l_level + 1
        else:
            if self.DIM == 1:
                l_level = np.array([np.linspace(0,self.level_max,self.level_max+1)]).transpose()
            else:
                l_level = allVL1leq(self.DIM, self.level_max)
        
        # filter out rows exceeding the individual level cap
        for i_DIM in range(self.DIM):
            l_level = l_level[l_level[:,i_DIM] <= self.level[i_DIM]]
        
        self.l_level = l_level
        
        # Consider interaction order (filter out multi-indices exceeding it)
        if self.interaction_order < self.DIM:
            if any("fejer2" in s for s in self.gridtype):
                self.l_level = self.l_level[np.sum(self.l_level>1,axis=1)<=self.interaction_order,:]
            else:
                self.l_level = self.l_level[np.sum(self.l_level>0,axis=1)<=self.interaction_order,:]
        
        # filter out inactive nodes multi-indices
        #self.l_level = self.l_level[np.sum(self.l_level,axis=1)==self.level_max,:]
        
        # make cubature lookup table for knots (dl_k) and weights (dl_w) [max(l) x DIM]
        print("    Generating difference grids ...")
        dl_k = [[0 for x in range(self.DIM)] for x in range(int(np.amax(level)+1))]
        dl_w = [[0 for x in range(self.DIM)] for x in range(int(np.amax(level)+1))]

        for i_dim in range(self.DIM):
    
            for i_level in level_sequence[i_dim]: #range(level[i_dim]+1):
 
                if self.gridtype[i_dim] == 'jacobi':    # Jacobi polynomials
                    knots_l, weights_l     = quadrature_jacobi_1D(order_sequence[i_dim][i_level],self.gridshape[0][i_dim]-1, self.gridshape[1][i_dim]-1)
                    knots_l_1, weights_l_1 = quadrature_jacobi_1D(order_sequence[i_dim][i_level-1],self.gridshape[0][i_dim]-1, self.gridshape[1][i_dim]-1)
                if self.gridtype[i_dim] == 'hermite':   # Hermite polynomials
                    knots_l, weights_l     = quadrature_hermite_1D(order_sequence[i_dim][i_level])
                    knots_l_1, weights_l_1 = quadrature_hermite_1D(order_sequence[i_dim][i_level-1])
                if self.gridtype[i_dim] == 'patterson': # Gauss-Patterson
                    knots_l, weights_l     = quadrature_patterson_1D(order_sequence[i_dim][i_level])
                    knots_l_1, weights_l_1 = quadrature_patterson_1D(order_sequence[i_dim][i_level-1])    
                if self.gridtype[i_dim] == 'cc':        # Clenshaw Curtis
                    knots_l, weights_l     = quadrature_cc_1D(order_sequence[i_dim][i_level])
                    knots_l_1, weights_l_1 = quadrature_cc_1D(order_sequence[i_dim][i_level-1])
                if self.gridtype[i_dim] == 'fejer2':    # Fejer type 2 
                    knots_l, weights_l     = quadrature_fejer2_1D(order_sequence[i_dim][i_level-1])
                    knots_l_1, weights_l_1 = quadrature_fejer2_1D(order_sequence[i_dim][i_level-2])       
                if i_level == 0 and (not self.gridtype[i_dim] == 'fejer2'):
                    dl_k[i_level][i_dim] = knots_l
                    dl_w[i_level][i_dim] = weights_l
                elif i_level == 1 and (self.gridtype[i_dim] == 'fejer2'): 
                    dl_k[i_level][i_dim] = knots_l
                    dl_w[i_level][i_dim] = weights_l
                else:
                    dl_k[i_level][i_dim] = np.hstack((knots_l, knots_l_1))
                    dl_w[i_level][i_dim] = np.hstack((weights_l, -weights_l_1))

        # make list of all tensor products according to multiindex list "l"
        print("    Generating subgrids ...")
        dL_k = []
        dL_w = []
        
        for i_l_level in range(self.l_level.shape[0]):
    
            knots_temp = []
            weights_temp = []
    
            for i_dim in range(self.DIM): 
        
                knots_temp.append(np.asarray(dl_k[np.int(self.l_level[i_l_level,i_dim])][i_dim],dtype=float))
                weights_temp.append(np.asarray(dl_w[np.int(self.l_level[i_l_level,i_dim])][i_dim],dtype=float))
    
            # tensor product of knots    
            dL_k.append(combvec(knots_temp))
    
            # tensor product of weights
            dL_w.append(np.prod(combvec(weights_temp),axis=1))

        dL_w = np.hstack(dL_w)    
        dL_k = np.vstack(dL_k)          

        # find similar points in grid and formulate Point list
        print("    Merging subgrids ...")
        Point_number_list = np.zeros(dL_w.shape[0])-1
        point_no = 0
        epsilon_k = 1E-6
        coords_norm= []

        while any(Point_number_list<0):
            notfound = Point_number_list < 0
            dL_k_nf = dL_k[notfound,:]
            point_temp = np.zeros(dL_k_nf.shape[0])-1
            point_temp[np.sum(np.abs(dL_k_nf-dL_k_nf[0,:]),axis=1)<epsilon_k] = point_no
            Point_number_list[notfound] = point_temp
            point_no = point_no + 1
            coords_norm.append(dL_k_nf[0,:])
        
        coords_norm = np.array(coords_norm)
        Point_number_list = np.asarray(Point_number_list,dtype=int)

        weights = np.zeros(np.amax(Point_number_list)+1)-999

        for i_point in range(np.amax(Point_number_list)+1):
            weights[i_point] = np.sum(dL_w[Point_number_list == i_point])

        # filter for very small weights
        print("    Filter grid for very small weights ...")
        epsilon_w = 1E-8/self.DIM
        keep_point = np.abs(weights) > epsilon_w
        self.weights = weights[keep_point]/2**self.DIM
        self.coords_norm = coords_norm[keep_point]
        self.level_sequence = level_sequence
        self.order_sequence = order_sequence
        
        # rescale normalized coordinates in case of normal distributions and "fejer2" or "cc" grids
        # +- 0.675 * sigma -> 50%
        # +- 1.645 * sigma -> 90%
        # +- 1.960 * sigma -> 95%        
        # +- 2.576 * sigma -> 99%
        # +- 3.000 * sigma -> 99.73%
        for i_DIM in range(self.DIM):
            if (self.pdftype[i_DIM] == "norm" or self.pdftype[i_DIM] == "normal") and (not(self.gridtype[i_DIM] == "hermite")): 
                self.coords_norm[:,i_DIM] = self.coords_norm[:,i_DIM]*1.960
                
        # denormalize grid to original parameter space
        print("    Denormalizing grid for computations ...")
        self.coords = denorm(self.coords_norm,self.pdftype,self.gridshape,self.limits)
        
#%% random-grid object subclass   
###############################################################################
  
class randomgrid():
    def __init__(self, pdftype, gridshape, limits, N):
        #print "-> Generate random grid"
        self.pdftype   = pdftype         # pdftype: "beta", "normal" [1 x DIM]
        self.gridshape = gridshape       # pdfshape: jacobi: -> [alpha and beta] hermite: -> [mean, variance] list [2 x DIM]
        self.limits    = limits          # limits: [min, max]  list [2 x DIM]  
        self.N         = int(N)          # Number of random samples
        self.DIM       = len(self.pdftype) # number of dimension
        
        #np.random.seed()  
        # generate random samples for each random input variable [N x DIM]
        self.coords_norm = np.zeros([self.N, self.DIM]) 
        
        for i_DIM in range(self.DIM):
            if self.pdftype[i_DIM] == "beta":
                self.coords_norm[:,i_DIM] = (np.random.beta(self.gridshape[0][i_DIM], self.gridshape[1][i_DIM], [self.N,1])*2.0 - 1)[:,0]
            if self.pdftype[i_DIM] == "norm" or self.pdftype[i_DIM] == "normal":
                self.coords_norm[:,i_DIM] = (np.random.normal(0, 1, [N,1]))[:,0]
                
        # denormalize grid to original parameter space
        #print "    Denormalizing grid for computations ..."
        self.coords = denorm(self.coords_norm,self.pdftype,self.gridshape,self.limits)

#==============================================================================
# #%% random-grid object subclass   
# ###############################################################################
# import pyDOE    
# class LHSgrid(grid):
#     def __init__(self, pdftype, gridshape, limits, N):
#         print "-> Generate LHS grid"
#         grid.__init__(self)
#         self.pdftype   = pdftype         # pdftype: "beta", "normal"
#         self.gridshape = gridshape       # pdfshape: jacobi: -> [alpha and beta] hermite: -> [mean, variance]
#         self.limits    = limits          # limits: [min, max]   
#         self.N         = int(N)          # Number of random samples
#         self.DIM       = len(self.pdftype) # number of dimension
#         
#         coords_raw = pyDOE.lhs(self.DIM, samples=N) # criterion='center'
#         self.coords_norm = np.zeros([self.N, self.DIM]) 
#         
#         for i_DIM in range(self.DIM):
#             if self.pdftype[i_DIM] == "beta":
#                 self.coords_norm[:, i_DIM] = scipy.stats.beta(self.gridshape[0][i_DIM], self.gridshape[1][i_DIM]).ppf(coords_raw[:, i_DIM])*2 - 1
#             if self.pdftype[i_DIM] == "norm" or self.pdftype[i_DIM] == "normal":
#                 self.coords_norm[:, i_DIM] = scipy.stats.norm(loc=0, scale=1).ppf(coords_raw[:, i_DIM])
#                 # denormalize grid to original parameter space
#         
#        # print "    Denormalizing grid for computations ..."
#         self.coords = denorm(self.coords_norm,self.pdftype,self.gridshape,self.limits)        
#==============================================================================
                
