#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 BinPo, TB-Poisson Solver for 2DES
 Copyright (C) 2021 BinPo Team
 
 BPModule.py is part of BinPo.
 
 BinPo is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 BinPo is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with BinPo. See ~/COPYING file or <https://www.gnu.org/licenses/>.
 
 DESCRIPTION:
    This is the BinPo main module.
"""

__author__ = "Emanuel A. Martinez"
__email__ = "emanuelm@ucm.es"
__copyright__ = "Copyright (C) 2021 BinPo Team"
__version__ = 1.1
__date__ = "August 9, 2022"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import interp1d
from scipy.special import expit
import scipy.linalg as LA
from ase.dft.kpoints import bandpath, monkhorst_pack
import numexpr as ne
from ase import Atoms
from ase.build import surface

# Physical constants
e0 = 8.854187817620389e-12 # vacuum permittivity in F/m = C/(V*m)
kB = 1.380649e-23 # Boltzmann constant in J/K
e = 1.602176634e-19 # elementary charge in C
eVtoJ = 1.602176634e-19 # Conversion from eV to J  

###############################################################################################################################
class PotentialEnergy:
   counter = 1

   def __init__(self, L, BC1 = None, BC2 = None, eps_model = None, V_init = 'linear'):
      """
      This class holds the potential energy along the slab in eV. Once an initial potential is chosen, 
      the same instance is updated at every step of Tight-Binding Poisson (TB-Poisson) scheme.

      Parameters
      ----------
      L : integer. Number of planes in the slab.
      BC1 : float, optional. Value of the potential energy in eV at the top-most layer.
      BC2 : float, optional. Value of the potential energy in eV at the bottom-most layer.
      eps_model : string, optional. Relative permittivity model.
      V_init : string, optional. Shape of the initial potential profile.
      """ 
      self.L = L
      self.BC1 = BC1
      self.BC2 = BC2
      self.eps_model = eps_model
      if V_init == 'linear' or V_init == 'lin':
          self.values = np.arange(self.L)*(self.BC2-self.BC1)/(self.L-1) + self.BC1 #start with triangular potential
      if V_init == 'exponential' or V_init == 'exp':
         decay = 2.0 # exponential decay constant
         self.values = self.BC1*np.exp(-np.arange(self.L)/decay) #start with exponential potential
         
#------------------------------------------------------------------------------------------------------------------------------
   def update_potential(self, new_values):
       """
       This method updates the values of potential energy.

       Parameters
       ----------
       new_values : array of float. New values of the potential at each plane.

       Returns
       -------
       PotentialEnergy instance with the values and counter updated.
       """
       self.values = new_values  
       self.counter += 1
       return self.values
   
#------------------------------------------------------------------------------------------------------------------------------
   def __repr__(self):
       """ Define a representation to print nicely the potential energy."""
       return 'V = ' + str(self.values)
   
#------------------------------------------------------------------------------------------------------------------------------
   def Error(self, Vnew):
       """
       Compute the error between potential energies.

       Parameters
       ----------
       Vnew : array of floats. Values of a new potential energy.

       Returns
       -------
       float. The error calculated between the new potential and the values of potential in the 
       current PotentialEnergy instance.
       """
       delta = (self.values - Vnew)**2/self.values[0]**2
       return np.sum(delta)/self.L
   
#------------------------------------------------------------------------------------------------------------------------------
   def to_tensor(self, tensor_size = 1, nwann = 6):
       """
       Method that transforms the potential energy array into an array of matrices. Each matrix in the 
       array has a dimension of 'nwann*L', while the dimension of the whole array is 'tensor_size'.

       Parameters
       ----------
       tensor_size : integer, optional. Dimension of the array of matrices.
       nwann : integer, optional. Number of elements in the MLWF basis.

       Returns
       -------
       ndarray of floats. An array of matrices containing the potential energy values ready to be added
       to a Hamiltonian. 
       """
       aux1 = np.array([self.values]*nwann)
       aux2 = np.ndarray.flatten(aux1.T)
       if tensor_size == 1:
           return np.identity(self.L*nwann)*aux2
       aux3 = np.identity(self.L*nwann)*aux2
       return np.array([aux3]*tensor_size)
   
#------------------------------------------------------------------------------------------------------------------------------   
   def electric_field(self, zdist, b = 2):
        """
        Compute the electric field from the potential energy.
 
        Parameters
        ----------
        zdist : float. Interplanar distance that discretizes the slab.
        b : integer, optional. Fraction of interplanar distance that defines the bin to compute the derivative.

        Returns
        -------
        array of floats. The electric field along the slab in V/m.
        """
        z0 = np.arange(self.L)
        # perform the gradient calculation
        f = interp1d(z0, self.values, kind = 'linear', fill_value = "extrapolate") # Function generated by linear interpolation
        z1 = np.arange(0, self.L, 1/b) # Auxiliary array to filter elements
        df = np.gradient(f(z1), 1/b, edge_order = 2) # gradient taking the first neighbors
        return df[np.arange(0, b*self.L, b)]/zdist
    
#------------------------------------------------------------------------------------------------------------------------------
   def permittivity_model(self, E):
       """
       Evaluate a specific model for the relative permittivity on the electric field profile
       along the slab, E(Rz). Predefined models can be chosen and any user-defined model can be applied
       by typing an expression with Python syntax from the configuration file scp.yaml.

       Parameters
       ----------
       E : array of floats. Electric field in V/m.

       Returns
       -------
       array of floats. Relative permittivity evaluated at each E(Rz)
       """
       if self.eps_model == 'Cop': # Copie model for STO(100) at 10 K
           X0 = 24000 
           Ec = 4.7e5
           return 1 + X0/(1 + E/Ec)

       elif self.eps_model == 'Ang': # Model: Fit of Ang data for KTO at 14 K
           X0 = 2837 
           Ec = 892244
           return 1 + X0/( 1 + (E/Ec)**2)**(2/5)
           
       elif 'cte' in self.eps_model:
           return np.ones(self.L)*float(self.eps_model.split('_')[1]) # Constant model
       else:
           try:
               eps = ne.evaluate(self.eps_model) # User-defined model
           except:
               raise ValueError('%s is an invalid function! Please, make sure of either selecting a permittivity from the available ones or correctly using Python syntax!' % self.eps_model)
           
           return eps
       
#------------------------------------------------------------------------------------------------------------------------------ 
   def Poisson_solver(self, rho, zdist, area, neumann_bc2 = False, background_Q = False, site_qdef = 0.0, ext_qdef = 5):
       '''
       Poisson equation solver. 

       Parameters
       ----------
       rho : array of float. Electron density in number of electrons/plane.
       zdist : float. Interplanar distance in m to define the 3D electron density.
       area : float. Area of 2D unit cell in m**2 to define the 3D electron density.
       neumann_bc2 : boolean, optional. Whether or not to use Neumann boundary condition at L-1.
       background_Q : boolean, optional. Whether or not to use a fixed background electron density.
       site_qdef : float, optional. Value of the fixed background electron density. Only used if 'background_Q = True'.
       ext_qdef : integer, optional. Extent of the fixed background electron density in planes. Only used if 'background_Q = True'.

       Returns
       -------
       Vout : array of float. Potential energy values in eV.
       n_iter : integer. The total number of iterations in Poisson step.
       '''
       err = 1 # Arbitrary initialization for the error between potentials
       n_iter = 0 # iteration index of Poisson equation
       Poisson_thr = 3e-12 # threshold for convergence in the Poisson equation
       Poisson_maxiter = 500000 # maximum number of iterations allowed for solving Poisson equation
       Vout = np.zeros_like(self.values)
       Vin = np.copy(self.values)
       
       rho_def = np.zeros(self.L)
       # The next block is only used if a fixed background density is set to True
       #------------------------------------------------------------------------------
       # creation of background of constant charge and extent = 'ext_qdef'
       if background_Q == True:
           qfix = np.ones(int(ext_qdef))*float(site_qdef)
           rho_def[:ext_qdef] = qfix
       #------------------------------------------------------------------------------ 
       while err > Poisson_thr: # the loop will run until convergence is achieved (or not)
           dV = self.electric_field(zdist)
           epz = self.permittivity_model(dV)
           
           Vout[0] = self.BC1 # setting Dirichlet boundary conditions
           if neumann_bc2 == False:
               Vout[self.L-1] = self.BC2
               
           for i in range(self.L):
               if i > 0 and i < self.L-1:
                   # interplanar volume model + constant background density
                   Vout[i] = 0.5*(Vin[i+1] + Vin[i-1] + e*(rho[i] + rho_def[i])*zdist/(area*e0*epz[i])) 
           
           if neumann_bc2 == True: # setting Neumann boundary condition if used
               Vout[self.L-1] = Vout[self.L-2]

           err = np.sum((np.array(Vout)-np.array(Vin))**2/self.BC1**2)/self.L # compute the error between Poisson iterations
           Vin = np.copy(Vout)
           n_iter += 1
           if n_iter == Poisson_maxiter: # number of allowed iterations in Poisson solver
               raise RuntimeError('Poisson step did not converge after ' + str(n_iter) + ' iterations.')
        
       return Vout, n_iter

###############################################################################################################################
class CrystalFeatures:

   def __init__(self, face, a0, material, c = 6.0, manifold = 't2g', sgeom = 'cubic'):
       """
       This class provides basic crystallographic characteristics of the system.

       Parameters
       ----------
       face : string. Confinement direction formatted as 'hkl'.
       a0 : float. Lattice parameter in Angs for the cubic or hexagonal unit-cell.
       c : float. Lattice parameter in Angs along c axis. Only needed for hexagonal unit-cell.
       material : string. Name of the material or system.
       sgeom : string. Unit-cell geometry.
       """
       if sgeom == 'cubic':
           if face in ['100','010','001']:
               self.face = '100'
           if face in ['110','101','011']:
               self.face = '110'
           if face == '111':
               self.face = '111'
       elif sgeom == 'hexagonal':
           if face in ['100','010','001']:
               self.face = '001h'
       else:
           raise ValueError("%s is an invalid system geometry. Allowed values are 'cubic' or 'hexagonal'." % sgeom)
       self.a0 = a0
       self.c = c
       self.material = material
       self.sgeom = sgeom
#------------------------------------------------------------------------------------------------------------------------------   
   def interplanar_distance(self):
       """ This method computes the interplanar distance in m."""
       if self.sgeom == 'cubic' and self.manifold == 't2g':
           a1, a2, a3 = np.fromiter(self.face, dtype = int)
           return self.a0/np.sqrt(a1*a1 + a2*a2 + a3*a3)*1e-10 # the factor 1e-10 is to convert Angs to m
       elif self.sgeom == 'cubic' and self.manifold == 'other':
           return self.a0*np.sqrt(3)*1e-10
       elif self.sgeom == 'hexagonal':
           return self.c*1e-10
#------------------------------------------------------------------------------------------------------------------------------
   def face_area(self):
       """ This method computes the unit cell area in m**2."""
       if self.sgeom == 'cubic':
           if self.face  in ['001','010','100']:
               return self.a0*self.a0*1e-20
           if self.face in ['110','101','011']:
               return np.sqrt(2)*self.a0*self.a0*1e-20 # the factor 1e-10 is to convert Angs**2 to m**2
           if self.face == '111':
               return 2*self.a0*self.a0*1e-20*(np.sqrt(3)/2)
       if self.sgeom == 'hexagonal':
           return self.a0*self.a0*(np.sqrt(3)/2)*1e-20
        
###############################################################################################################################
class Quasi2DHamiltonian():   
    # setting arbitrarily some general attributes
    T = 10 # Temperature in K
    ef = 11.5997 # Fermi level
    kpts = 15*15 # Total numbers of k-points
    HOL = -np.inf # Highest occupied level
    
    def __init__(self, dix, L, nwann):
        """
        This class represents the slab Hamiltonian with the operations involved in 
        the TB step of the algorithm.        
        
        Parameters
        ----------
        dix: dictionary containing the nwann*nwann complex matrices obtained from
        the WF separation for each plane.
        L : integer. Number of planes in the slab.
        nwann : integer. Number of elements in the MLWF basis.
        """
        self.dix = dix
        self.L = L
        self.nwann = nwann
#------------------------------------------------------------------------------------------------------------------------------
    @classmethod
    def set_parameters(cls, Temp, Ef, Kpts, HOL):
        """ Set the values of the general attributes for this class."""
        cls.T = Temp
        cls.ef = Ef
        cls.kpts = Kpts
        cls.HOL = HOL
#------------------------------------------------------------------------------------------------------------------------------
    def HamiltonianTensor(self, nwann):
        """
        This method creates the slab Hamiltonian matrices for each k-point in the grid and arranges it in an array.
        Parameters
        ----------
        nwann : integer. Number of elements in the MLWF basis.
        Returns
        -------
        ndarray of complex. Hamiltonian tensor, namely, an array of lenght Nk*Nk containing the slab Hamiltonian matrices,
        with 'Nk' being the square root of the total k-points.
        """
        GL = []
        for j in range(self.kpts): # loop over total k-points
            Z = np.reshape(self.dix['T_Z-7'][j],(nwann, nwann)) #---------> z = -7 (interactions between planes -7 and 0)
            A = np.reshape(self.dix['T_Z-6'][j],(nwann, nwann))
            B = np.reshape(self.dix['T_Z-5'][j],(nwann, nwann))
            C = np.reshape(self.dix['T_Z-4'][j],(nwann, nwann))
            D = np.reshape(self.dix['T_Z-3'][j],(nwann, nwann))
            E = np.reshape(self.dix['T_Z-2'][j],(nwann, nwann))
            F = np.reshape(self.dix['T_Z-1'][j],(nwann, nwann))
            G = np.reshape(self.dix['T_Z0'][j],(nwann, nwann)) #---------> z = 0 (interactions within the plane 0)
    
            ZZ = np.kron(np.eye(self.L,k=-7),Z) #---------> z = -7
            AA = np.kron(np.eye(self.L,k=-6),A)
            BB = np.kron(np.eye(self.L,k=-5),B)
            CC = np.kron(np.eye(self.L,k=-4),C)
            DD = np.kron(np.eye(self.L,k=-3),D)
            EE = np.kron(np.eye(self.L,k=-2),E)
            FF = np.kron(np.eye(self.L,k=-1),F)
            GG = np.kron(np.eye(self.L,k=0),G) #---------> z = 0
        
            GL.append(ZZ + AA + BB + CC + DD + EE + FF + GG) # Append the total Hamiltonian matrix at each k-point 
        return np.array(GL)
#------------------------------------------------------------------------------------------------------------------------------   
    @staticmethod
    def SumInCell(X,l):
        """ Sum the square modulus of the eigenvectors per site."""
        return np.linalg.norm(np.array(np.split(X,l)),axis = 1)**2
#------------------------------------------------------------------------------------------------------------------------------        
    @classmethod
    def FermiOcc(cls, E):
        """ Fermi-Dirac occupation function."""
        return expit(-1*eVtoJ*(E-cls.ef)/(kB*cls.T))    
#------------------------------------------------------------------------------------------------------------------------------
    @classmethod
    def ChargeSolver(cls, Hk, V, nwann, kBT_ext = 3):
        """
        This method is the charge equation solver used in case. It diagonalizes the Hamiltonian matrices
        and computes the electron density from the eigenvectors. 

        Parameters
        ----------
        Hk : ndarray of complex. The Hamiltonian tensor.
        V : A PotentialEnergy instance.
        nwann : integer. Number of elements in the MLWF basis.
        kBT_ext : integer, optional. Number of how many k_B*T factors consider to set the upper limit
        in energy for diagonalization and charge density calculation.

        Returns
        -------
        rho : array of floats. The electron density in number of electrons/plane.
        """
        Hv = V.to_tensor(nwann = nwann)
        deltaE = kBT_ext*eVtoJ*kB*cls.T # writting the kBT_ext in eV
        fLim = 1.01 # this factor increases a bit the value of HOL to set a reference
                    # lower limit in energy for diagonalization and charge density calculation.
        Lrho = []
        for hk in Hk: 
            HT = hk + Hv
            # It tries to use scipy routine, which is faster. Otherwise, numpy routine will be used.
            try:
                # with subset_by_value we filter the range of interest, which is [HOL*fLim, ef + deltaE]
                w, v = LA.eigh(HT, lower = True, overwrite_a = True, subset_by_value = (cls.HOL*fLim, cls.ef + deltaE))
            except:
                w0, v = np.linalg.eigh(HT, UPLO = 'L')
                w = w0[np.where(np.logical_and(w0 > cls.HOL, w0 < (cls.ef + deltaE)))]
            
            # Append the value of the electron density along the slab for a specific k-point
            Lrho.append(np.dot(cls.FermiOcc(w).T, np.array([cls.SumInCell(v[:, j],V.L) for j in range(len(w))])))
        # Finally, values are averaged over all k-points in grid    
        rho = np.sum(np.array(Lrho, dtype = object), axis = 0)/cls.kpts # Charge density in e/plane
        Lrho.clear()

        return rho
#------------------------------------------------------------------------------------------------------------------------------    
    def HamiltonianMatrix(self, i):
        """
        This is an auxiliary method for ChargeSolver2, which is applied if Total_Hk_method = 'iterable'.
        It generates the Hamiltonian matrix for just one k-point at a time.

        Parameters
        ----------
        i : integer. Index of the k-point to consider.

        Returns
        -------
        ndarray of complex. Hamiltonian matrix in a specific k-point.

        """
        Z = np.reshape(self.dix['T_Z-7'][i],(self.nwann, self.nwann)) #---------> z = -7
        A = np.reshape(self.dix['T_Z-6'][i],(self.nwann, self.nwann))
        B = np.reshape(self.dix['T_Z-5'][i],(self.nwann, self.nwann))
        C = np.reshape(self.dix['T_Z-4'][i],(self.nwann, self.nwann))
        D = np.reshape(self.dix['T_Z-3'][i],(self.nwann, self.nwann))
        E = np.reshape(self.dix['T_Z-2'][i],(self.nwann, self.nwann))
        F = np.reshape(self.dix['T_Z-1'][i],(self.nwann, self.nwann))
        G = np.reshape(self.dix['T_Z0'][i],(self.nwann, self.nwann)) #---------> z = 0
    
        ZZ = np.kron(np.eye(self.L,k=-7),Z) #---------> z = -7
        AA = np.kron(np.eye(self.L,k=-6),A)
        BB = np.kron(np.eye(self.L,k=-5),B)
        CC = np.kron(np.eye(self.L,k=-4),C)
        DD = np.kron(np.eye(self.L,k=-3),D)
        EE = np.kron(np.eye(self.L,k=-2),E)
        FF = np.kron(np.eye(self.L,k=-1),F)
        GG = np.kron(np.eye(self.L,k=0),G) #---------> z = 0
        
        return np.array(ZZ + AA + BB + CC + DD + EE + FF + GG)
#------------------------------------------------------------------------------------------------------------------------------   
    @classmethod
    def ChargeSolver2(cls, H, V, n_wann, kBT_ext = 3):
        '''
        Alternative and less memory-consuming method to get the electron density.
        Used when Total_Hk_method = 'iterable'.

        Parameters
        ----------
        H : A Quasi2DHamiltonian instance. 
        V : A PotentialEnergy instance.
        nwann : integer. Number of elements in the MLWF basis.
        kBT_ext : integer, optional. Number of how many k_B*T consider to set the upper limit
        in energy for diagonalization and charge density calculation.

        Returns
        -------
        rho : array of floats. The electron density in number of electrons/plane.
        '''
        Hv = V.to_tensor(nwann = n_wann)
        deltaE = kBT_ext*eVtoJ*kB*cls.T # writting the kBT_ext in eV
        fLim = 1.01 # this factor increases a bit the value of HOL to set a reference
                    # lower limit in energy for diagonalization and charge density calculation.
        Lrho = []
        for i in range(cls.kpts): # loop over the k-grid
            HT = H.HamiltonianMatrix(i) + Hv # construction of Hamiltonian matrix
            # It tries to use scipy routine, which is faster. Otherwise, numpy routine will be used.
            try:
                # with subset_by_value we filter the range of interest, which is [HOL*fLim, ef + deltaE]
                w, v = LA.eigh(HT, lower = True, overwrite_a = True, subset_by_value = (cls.HOL*fLim, cls.ef + deltaE))
            except:
                w0, v = np.linalg.eigh(HT, UPLO = 'L')
                w = w0[np.where(np.logical_and(w0 > cls.HOL, w0 < (cls.ef + deltaE)))]
            
            # Append the value of the electron density along the slab for specific k-point
            Lrho.append(np.dot(cls.FermiOcc(w).T, np.array([cls.SumInCell(v[:, j],V.L) for j in range(len(w))])))
        # Finally, values are averaged over all k-points in grid    
        rho = np.sum(np.array(Lrho, dtype = object), axis = 0)/cls.kpts # Charge density in e/plane
        Lrho.clear()
        return rho
    
###############################################################################################################################

def Wann_Sep(Zn, nwann):
    """
    This function separate the r-space bulk Hamiltonian for a particular plane z into the expectation values 
    with the MLWF indices w, w'; i. e., it separates <0w|H|Rxy + Rz w'> matrix elements.

    Parameters
    ----------
    Zn : string. Name of file containing the r-space bulk Hamiltonian in a particular plane.
    nwann : integer. Number of elements in the MLWF basis.
    Returns
    -------
    D : dictionary of {(w, w'), <0w|H|Rxy + Rz w'>} pairs.
    """
    D = {}
    for i in range(nwann):
        for j in range(nwann):  
            D[str(i+1) + '-' + str(j+1)] = []
    
    for z in Zn:
        D[str(int(z[3])) + '-' + str(int(z[4]))].append(z)

    return D
#------------------------------------------------------------------------------------------------------------------------------
def Hopping2D(kmesh, D):
    """
    This function executes the 2D Fourier transform for a particular plane Rz.

    Parameters
    ----------
    kmesh : ndarray. 2D k-grid.
    D : dictionary of {(w, w'), <0w|H|Rxy + Rz w'>} pairs in real space, with w,w' being the MLWF indices.

    Returns
    -------
    ndarray of complex. Array with the k-space transformed <0w|H|Rxy + Rz w'> elements.
    """
    d = D.copy() # Create copy of input dictionary
    kx = kmesh.T[0] # Take the x and y grid components in separated arrays
    ky = kmesh.T[1]
    
    for key,val in d.items(): # Here, the 2D Fourier transform is computed
        tK = np.zeros_like(kx, dtype = np.complex64)
        for v in val:
            tK += (v[5] + 1j*v[6])*np.exp(-1j*2*np.pi*(kx*v[0] + ky*v[1]))
        d[key] = tK # The k-transformed value is associated to the new dictionary
    
    # Create an array from the new dictionary 
    L = []
    for val in d.values():
        L.append(np.array(val))
         
    SM = np.split(np.array(L).T, len(kx)) # Split the new array to match the numbers of k-points
    
    return np.array(SM)

#------------------------------------------------------------------------------------------------------------------------------
def Kmeshgrid(Nk = 15, scale = 1.0, delta_kx = 0.0, delta_ky = 0.0):
    """
    It generates the k-grid using the Monkhorst-Pack method as implemented in ASE. A shift and a 
    scale factor to the grid can be set too.

    Parameters
    ----------
    Nk : integer, optional. Square root of the total k-points in the grid.
    scale : float, optional. Factor to affect the area encompassed by the grid.
    delta_kx : float, optional. Shift in k_x direction of the grid.
    delta_ky : float, optional. Shift in k_y direction of the grid.

    Returns
    -------
    ndarray of floats. Array containing every k-point of the grid. 
    """
    delta = np.array((delta_kx, delta_ky, 0.0)) # define the shift in k_x, k_y
    return monkhorst_pack((Nk,Nk,1))*scale + delta

#------------------------------------------------------------------------------------------------------------------------------
def Kmeshgrid2(Nk = 15, scale = 1.0, delta_kx = 0.0, delta_ky = 0.0, a = 1.0):
    """
    It generates the k-grid using the Monkhorst-Pack method as implemented in ASE. A shift and a 
    scale factor to the grid can be set too. It is a variation of Kmeshgrid to be used in energy
    slices calculations.

    Parameters
    ----------
    Nk : integer, optional. Square root of the total k-points in the grid.
    scale : float, optional. Factor to affect the area encompassed by the grid.
    delta_kx : float, optional. Shift in k_x direction of the grid.
    delta_ky : float, optional. Shift in k_y direction of the grid.
    a = float, optional. Lattice parameter in Angs.
    Returns
    -------
    ndarray of floats. Array containing every k-point of the grid. 
    """
    delta = np.array((delta_kx, delta_ky, 0.0)) # define the shift in k_x, k_y
    return monkhorst_pack((Nk,Nk,1))*scale + delta*a/(2*np.pi)

#------------------------------------------------------------------------------------------------------------------------------
def BandPath(path_str, cell, kpoints):
    """
    It generates a path within the irreducible first Brillouin zone (IBZ1).

    Parameters
    ----------
    path_str : string. High-symmetry points in IBZ1 to draw the path.
    cell : ASE Cell object. The 2D unit cell involved.
    kpoints : integer. Number of k-points to discretize the path.

    Returns
    -------
    ASE BandPath instance. The bandpath taken in IBZ1.
    """
    return bandpath(path_str, cell, npoints = kpoints)
 
#------------------------------------------------------------------------------------------------------------------------------
def CellGenerator(material, face, a0):
    """
    Function that generates the 2D unit cell for a system depending on the crystal face and geometry.

    Parameters
    ----------
    material : string. Name of material or system. 
    face : string. The lower indices 'hkl' crystal face where the surface is created.
    a0 : float. Lattice parameter as obtained from ab initio calculations.

    Returns
    -------
    cell : ASE Cell object. The 2D unit cell involved.
    """
    if material == 'STO' or material == 'STOB' or material == 'STOC':
        # created the cubic STO perovskite
        STO = Atoms('SrTiO3',
                cell = [a0,a0,a0],
                pbc = True,
                scaled_positions = ((0,0,0),(0.5,0.5,0.5),(0.5,0.5,0),(0.5,0,0.5),(0,0.5,0.5))
                )
        # create a Cell ASE object according to the orientation 'hkl'
        if face in ['100','010','001']:
            STO001 = surface(STO, (0, 0, 1), 1, periodic = False)
            return STO001.cell
        if face in ['110','101','011']:
            STO110 = surface(STO, (1, 1, 0), 2, periodic = False)
            return STO110.cell            
        if face == '111':
            STO111 = surface(STO, (1, 1, 1), 3, periodic = False)
            return STO111.cell        
    
    elif material == 'KTO':
        # created the cubic KTO perovskite
        KTO = Atoms('KTaO3',
                cell = [a0,a0,a0],
                pbc = True,
                scaled_positions = ((0,0,0),(0.5,0.5,0.5),(0.5,0.5,0),(0.5,0,0.5),(0,0.5,0.5))
                )
        # create a Cell ASE object according to the orientation 'hkl'
        if face in ['100','010','001']:
            KTO001 = surface(KTO, (0, 0, 1), 1, periodic = False)
            return KTO001.cell
        if face in ['110','101','011']:
            KTO110 = surface(KTO, (1, 1, 0), 2, periodic = False)
            return KTO110.cell            
        if face == '111':
            KTO111 = surface(KTO, (1, 1, 1), 3, periodic = False)
            return KTO111.cell
            
    elif material == 'BTB':
        BTB = Atoms('BiTeBr',
                    cell = [4.2662, 4.2662, 6.487, 90.0, 90.0, 120.0],
                    pbc = True,
                    scaled_positions=((0.33, 0.67, 0.93), (0.00, 0.00, 0.658), (0.67, 0.33, 0.230))
                    )
        if face == '001h':
            BTB001 = surface(BTB, (0,0,1), 1, periodic = False)
            return BTB001.cell

    # Add a new material here if another W90 file is available!
    # Just copy an edit one of the above blocks

#------------------------------------------------------------------------------------------------------------------------------        
def SnapPlotter(fig, ax, V, Vin, Vout, delta_t, err, pause):
    """
    Function to generate a step by step visualization of the potential profile evolution.

    Parameters
    ----------
    fig, ax : Matplotlib Figure and Axes instances.
    V : Potential instance. Current Potential instance to get information for the plot.
    Vin : array of floats. Values of the potential input in TB-Poisson scheme at a specific iteration.
    Vout : array of floats. Values of the potential output in TB-Poisson scheme at a specific iteration.
    delta_t : float. Time spent in seconds.
    err : float. Total error in the TB-Poisson scheme.
    pause : float. Waiting time in seconds while showing the plot.
    """
    plt.ion()
    ax.clear()
    ax.plot(np.arange(V.L), Vin, marker = 'o', label = 'V$_{in}$')
    ax.plot(np.arange(V.L), Vout, marker = 'o', label='V$_{out}$')
    ax.set_title('ITERATION #' + str(V.counter) + '\nTotal error = ' + "{:.8f}".format(err) + ", Time spent = " + "{:.2f}".format(delta_t) + ' seg', size = 12)
    ax.hlines(0,-5,V.L+5, ls = ':', color = 'k')
    ax.set_xlabel('planes', size = 14)
    ax.set_ylabel('V(z) [eV]', size = 14)
    ax.set_xlim(-1,V.L)
    #ax.set_ylim(V.BC1*1.2,abs(V.BC1)/2) # uncomment this line to let the y-axis be fixed during the iterations
    ax.legend(loc = 'lower right', facecolor = 'w', framealpha = 0.6, fontsize = 'x-large')
    fig.canvas.draw()
    plt.pause(pause)
 
#------------------------------------------------------------------------------------------------------------------------------
def SnapError(fig, ax, lerr, tol, pause):
    """
    Function to generate a step by step visualization of the total error as function of the iterations.

    Parameters
    ----------
    fig, ax : Matplotlib Figure and Axes instances.
    lerr : array of floats. Array containing the total error values for each iteration.
    tol : float. Convergence threshold for the self-consistent solution.
    pause : float. Waiting time in seconds while showing the plot.
    """
    plt.ion()
    ax.clear()    
    ax.plot(np.arange(len(lerr)) + 1, np.array(lerr), marker = 's', color = 'k')
    ax.set_xlim(0, len(lerr)+1)
    ax.set_yscale('log')
    ax.set_title('Total error vs iterations', size = 10)
    ax.hlines(tol,-2,500, ls = 'dashed', color = 'r')
    fig.canvas.draw()
    plt.pause(pause)

#------------------------------------------------------------------------------------------------------------------------------
def orb_char(eigvec, orb, L, nwann = 6):
    """
    Function to compute the orbital character of the eigenvectors.
    At present, this is only implemented with cubic system and t2g manifold,
    like STO and KTO.

    Parameters
    ----------
    eigvec : array of complex. Eigenvector where the orbital character will be computed.
    orb : string. Orbitals to be considered.
    L : integer. Number of planes in the slab.
    nwann : integer, optional. Number of elements in the MLWF basis.

    Returns
    -------
    float. The fraction of the orbital character of the eigenvector.
    """
    if orb == 'zx':
        zx_u = np.arange(0,nwann*L,nwann) # The indices chosen to extract each orbital character depends on the
        zx_d = np.arange(1,nwann*L,nwann)  # original Wannier interpolation.
        zx = np.sort(np.concatenate([zx_u,zx_d]))
        eig_orb = eigvec[zx]
        
    if orb == 'yz':
        yz_u = np.arange(2,nwann*L,nwann)
        yz_d = np.arange(3,nwann*L,nwann)
        yz = np.sort(np.concatenate([yz_u,yz_d]))
        eig_orb = eigvec[yz]
      
    if orb == 'xy':
        xy_u = np.arange(4,nwann*L,nwann)
        xy_d = np.arange(5,nwann*L,nwann)
        xy = np.sort(np.concatenate([xy_u,xy_d]))
        eig_orb = eigvec[xy]
       
    return np.linalg.norm(eig_orb)**2
#-------------------------------------------------------------------------------------------------------------------------------
def BandCalculation_t2g(HK, HV, nbands = 50): 
    """
    Compute the total bandstructure for a specific slab Hamiltonian matrix.
    For cubic t2g manifold systems like STO and KTO.

    Parameters
    ----------
    HK : ndarray of complex. Hamiltonian tensor for a discretized path in IBZ1.
    HV : ndarray of floats. Potential energy tensor of equal dimension of HK.
    nbands : integer, optional. Number of bands to consider in the calculation.

    Returns
    -------
    ndarray of floats. A 2D array containing all the corresponding eigenvalues.
    """
    Ek = []
    for hk in HK:
        HT = hk + HV
        w, v = np.linalg.eigh(HT, UPLO='L')
        Ek.append(np.real(w[:nbands]))
    return np.array(Ek)
#-------------------------------------------------------------------------------------------------------------------------------
def BandCalculation_gen(HK, HV, skipbands = 0, nbands = 50): 
    """
    Compute the total bandstructure for a specific slab Hamiltonian matrix.
    For cubic or hexagonal systems and arbitrary manifolds.

    Parameters
    ----------
    HK : ndarray of complex. Hamiltonian tensor for a discretized path in IBZ1.
    HV : ndarray of floats. Potential energy tensor of equal dimension of HK.
    skipbands: integer, optional. Number of lower energy bands to skip.
    nbands : integer, optional. Number of bands to consider in the calculation.

    Returns
    -------
    ndarray of floats. A 2D array containing all the corresponding eigenvalues.
    """
    Ek = []
    for hk in HK:
        HT = hk + HV
        w, v = np.linalg.eigh(HT, UPLO='L')
        Ek.append(np.real(w[skipbands:skipbands + nbands]))

    return np.array(Ek)
#--------------------------------------------------------------------------------------------------------------------------------
def BandCalc_OrbitalChar_t2g(HK, HV, L, nbands = 50): 
    """
    Compute the orbital projected bandstructure for a specific slab Hamiltonian matrix.
    For cubic t2g manifold systems like STO and KTO.

    Parameters
    ----------
    HK : ndarray of complex. Hamiltonian tensor for a discretized path in IBZ1.
    HV : ndarray of floats. Potential energy tensor of equal dimension of HK.
    L : integer. Number of planes in the slab.
    nbands : integer, optional. Number of bands to consider in the calculation.

    Returns
    -------
    tuple of ndarrays of floats. Return the eigenvalues with the three components of the orbital character.
    """
    Ek = []
    Oc_zx, Oc_yz, Oc_xy = [], [], []
    OCzx, OCyz, OCxy = [], [], []
    for hk in HK:
        HT = hk + HV
        w, v = np.linalg.eigh(HT, UPLO='L')
        Ek.append(np.real(w[:nbands]))       
        for i in range(nbands):
            Oc_zx.append(orb_char(v[:, i], 'zx', L))
            Oc_yz.append(orb_char(v[:, i], 'yz', L))
            Oc_xy.append(orb_char(v[:, i], 'xy', L))
        OCzx.append(np.array(Oc_zx))
        OCyz.append(np.array(Oc_yz))
        OCxy.append(np.array(Oc_xy))
        Oc_zx.clear()
        Oc_yz.clear()
        Oc_xy.clear()
    return np.array(Ek), np.array(OCxy), np.array(OCyz), np.array(OCzx) 

#------------------------------------------------------------------------------------------------------------------------------
def PlaneProjector_t2g(HK, HV, L, i_plane, f_plane, nbands = 50):
    """
    Compute the planes projected bandstructure for a specific slab Hamiltonian matrix.
    For cubic t2g manifold systems like STO and KTO.

    Parameters
    ----------
    HK : ndarray of complex. Hamiltonian tensor for a discretized path in IBZ1.
    HV : ndarray of floats. Potential energy tensor of equal dimension of HK.
    L : integer. Number of planes in the slab.
    i_plane : integer. Plane index from which starts the projection.
    f_plane : integer. Plane index to stop the projection.
    nbands : integer, optional. Number of bands to consider in the calculation.

    Returns
    -------
    tuple of ndarrays of floats. Return the eigenvalues and the planes projection values.
    """
    Ek, D_oC, Oc = [], [], []
    for hk in HK:
        HT = hk + HV
        w, v = np.linalg.eigh(HT, UPLO='L')
        Ek.append(np.real(w[:nbands]))
        for i in range(nbands):
           Oc.append(np.sum(Quasi2DHamiltonian.SumInCell(v[:, i], L)[i_plane:f_plane]))
        D_oC.append(np.array(Oc))
        Oc.clear()
    return np.array(Ek), np.array(D_oC)
#------------------------------------------------------------------------------------------------------------------------------
def PlaneProjector_gen(HK, HV, L, i_plane, f_plane, skipbands = 0, nbands = 50):
    """
    Compute the planes projected bandstructure for a specific slab Hamiltonian matrix.
    For cubic or hexagonal systems and arbitrary manifolds.

    Parameters
    ----------
    HK : ndarray of complex. Hamiltonian tensor for a discretized path in IBZ1.
    HV : ndarray of floats. Potential energy tensor of equal dimension of HK.
    L : integer. Number of planes in the slab.
    i_plane : integer. Plane index from which starts the projection.
    f_plane : integer. Plane index to stop the projection.
    skipbands: integer, optional. Number of lower energy bands to skip.
    nbands : integer, optional. Number of bands to consider in the calculation.

    Returns
    -------
    tuple of ndarrays of floats. Return the eigenvalues and the planes projection values.
    """
    Ek, D_oC, Oc = [], [], []
    for hk in HK:
        HT = hk + HV
        w, v = np.linalg.eigh(HT, UPLO='L')
        Ek.append(np.real(w[skipbands: skipbands + nbands]))
        for i in range(nbands):
           Oc.append(np.sum(Quasi2DHamiltonian.SumInCell(v[:, i], L)[i_plane:f_plane]))
        D_oC.append(np.array(Oc))
        Oc.clear()
    return np.array(Ek), np.array(D_oC)

#------------------------------------------------------------------------------------------------------------------------------
def Maxwell_triangle(color1, color2, color3, fig_ax, rel_size = "30%", location = 4, border = 1, fontsize = 15):  
    """
    Function for generating a Maxwell color triangle. For the time being, used in
    cubic t2g manifold systems like STO and KTO.

    Parameters
    ----------
    color1, color2, color3 : arrays with the values of R, G, B components.
    fig_ax : Matplotlib Axes instance. Axis where the triangle is plotted.
    rel_size : string, optional. Percentage of the 'fig_ax' axis that is used by the triangle. 
    location : integer, optional. Where the triangle is placed.
    border : float, optional. Relative distance to the border of 'fig_ax' axis.
    fontsize : integer, optional. Size of the text placed at the vertices of the triangle.
    """
    ax_in = inset_axes(fig_ax, width = rel_size, height = rel_size, loc = location, borderpad = border)
    Nlines, Ncol = 150, 150
    img = np.zeros((Nlines,Ncol,4))
    dx, dy = 2.0/(Ncol-1), 1.0/(Nlines-1)
    for i in range(Ncol-1):
        for j in range(Nlines-1):
            x, y = -1.0+i*dx, j*dy
            g = y
            r = (x+1-g)/2.0
            b = 1.0-g-r
            if (r>=0) and (r<=1.0) and (g>=0) and (g<=1.0) and (b>=0) and (b<=1.0):
                r_c1 = r*cm.to_rgba_array(color1)
                g_c2 = g*cm.to_rgba_array(color2)
                b_c3 = b*cm.to_rgba_array(color3)
                aux = np.ndarray.flatten(r_c1 + g_c2 + b_c3)
                add_C = aux/np.max(aux[:3])
                add_C[3] = 1.0
                img[j][i] = add_C 
            else:
                img[j][i] = np.array([1.0, 1.0, 1.0, 0.0])
    a = 1.0/np.sqrt(3)
    ax_in.plot([a, 0, -a, a],[0, 1, 0, 0], 'w-', lw = 1.5)
    ax_in.imshow(img, origin = 'lower', interpolation = 'bilinear', extent = [-a, a, 0.0, 1.0])
    ax_in.annotate('xy', (a+0.02, -0.15), c = color1, size = fontsize)
    ax_in.annotate('yz', (-0.09, 1.05), c = color2, size = fontsize)
    ax_in.annotate('zx', (-a-0.35, -0.15), c = color3, size = fontsize)
    plt.axis([-1, 1, -0.2, 1.2])
    plt.axis('off')  
    return None

#------------------------------------------------------------------------------------------------------------------------------
def Gen_Colors(colors_rgb, color1, color2, color3, Norm = True):
    """
    Function to generate a linear combination of colors through a R,G,B trio.

    Parameters
    ----------
    colors_rgb : ndarray of floats. 
    color1, color2, color3 : strings. Matplotlib color name that will become the new color basis.
    Norm : boolean, optional. Whether or not to normalize the color trio values.

    Returns
    -------
    ndarray of floats. The colors_rgb ndarray transformed into the new color basis.
    """
    Lcol, Ncol = [], []
    for c in colors_rgb:
        Lcol.append(np.array(c[0]*np.array(cm.to_rgb(color1))+c[1]*np.array(cm.to_rgb(color2))+c[2]*np.array(cm.to_rgb(color3))))
    if Norm == True:
        for i in range(len(Lcol)):
            #Ncol.append(Lcol[i]/np.linalg.norm(Lcol[i])) # uncomment for not applying contrast enhancement
            Ncol.append(Lcol[i]/np.max(Lcol[i])) # normalization with contrast enhancement
        return np.array(Ncol)
    else:
        return np.array(Lcol) 

###############################################################################################################################

if __name__ == "__main__":
    print('\n')
    print('=========================================================================')
    print('=========================================================================')
    print('                                BinPo                                    ')
    print('=========================================================================')
    print('=========================================================================')
    print('  WELCOME TO THE TIGHT-BINDING POISSON CODE FOR 2DES!')
    print('\n')
    print('Author: Emanuel A. Martinez, Universidad Complutense de Madrid, Spain.')
    print('\n')
    print('BinPo main module. This is not a runnable file.')
    print('\n')
    print('=========================================================================')
    print('\t DOCUMENTATION')
    print('=========================================================================')
    print(help(__name__))    
    print('=========================================================================')
    print('=========================================================================')
    
