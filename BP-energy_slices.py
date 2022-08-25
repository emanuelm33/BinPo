#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 BinPo, TB-Poisson Solver for 2DES
 Copyright (C) 2021 BinPo Team
 
 BP-energy_slices.py is part of BinPo.
 
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
     ENERGY SLICES CALCULATION
     BP-energy_slices.py computes a constant energy plot of the 
     band structure.
"""

__author__ = "Emanuel A. Martinez"
__email__ = "emanuelm@ucm.es"
__copyright__ = "Copyright (C) 2021 BinPo Team"
__version__ = 1.1
__date__ = "August 9, 2022"

import numpy as np
import BPmodule as BPM
import time as t
import datetime as dt
import argparse
import yaml
import scipy.linalg as LA
import logging
import os

# Loading the configuration files to set parameters not defined by terminal
#--------------------------------------------------------------------------------------------------------------------------
with open('./config_files/energy_slices.yaml', 'r') as f:
    try:
        data = yaml.load(f, Loader=yaml.FullLoader)
    except yaml.YAMLError as exc:
        print ("Error in configuration file: ", exc)
#--------------------------------------------------------------------------------------------------------------------------
# Messages to print by terminal when checking the description with -h or --help command
msj_f = 'Identifier to compute the energy slice. String.'
msj_ec = 'Energy at which the slice is computed. Float. 0.0 is the Fermi level.'
msj_nk = 'Square root of the total number of points in 2D k-grid. Integer.'
msj_ba = 'Number of batches to split the total k-points. Total k-points per batch will be '\
    'sqrt_kgrid_numbers**2/batches. In consequence, batches has to be divisor '\
        'of sqrt_kgrid_numbers**2, otherwise an error will arise.'
msj_bf = 'Factor to modify the box limits in k-space. Float. Value 1.0 corresponds to the BZ1.'
msj_dk = 'K-grid offset. Float (tuple). (0.0, 0.0) corresponds to null offset.'
description = "ENERGY SLICES CALCULATION"
epilog = "By default, arguments not specified are taken from './config_files/energy_slices.yaml'."
#--------------------------------------------------------------------------------------------------------------------------
# Passing arguments by terminal and setting their default values
parser = argparse.ArgumentParser(description = description, epilog = epilog)
parser.add_argument('-id', type = str, default = data['ENERGY_SLICES']['identifier'], help = msj_f)
parser.add_argument('-ec', type = float, default = data['ENERGY_SLICES']['energy_cut'], help = msj_ec)
parser.add_argument('-nk', type = int, default = data['ENERGY_SLICES']['sqrt_kgrid_numbers'], help = msj_nk)
parser.add_argument('-ba', type = int, default = data['ENERGY_SLICES']['batches'], help = msj_ba)
parser.add_argument('-bf', type = float, default = data['ENERGY_SLICES']['kbox_factor'], help = msj_bf)
parser.add_argument('-dk', nargs = 2, type = float, default = data['ENERGY_SLICES']['kbox_shift'], help = msj_dk)
args = parser.parse_args()
#--------------------------------------------------------------------------------------------------------------------------
# Setting at once what depends on the identifier
identifier = args.id # Identifier for the calculation
# Loading the .YAML file created after SCP calculation
with open(identifier + '/' + identifier + '.yaml', 'r') as f:
    try:
        params = yaml.load(f, Loader=yaml.FullLoader)
    except yaml.YAMLError as exc:
        print ("Error in configuration file: ", exc)
# Loading the SC solution file
TBP = np.loadtxt(identifier + '/' + identifier + '_SCP.dat')
#--------------------------------------------------------------------------------------------------------------------------
# Copying the parameters updatable by terminal and loading
# the remaining parameters from energy_slices.yaml
Ecut = args.ec
Nk = args.nk
Nfr = args.ba
kbox = args.bf
dk_x, dk_y = args.dk
wec = data['ENERGY_SLICES']['win_energy_calc']
outfile = data['ENERGY_SLICES']['outfile']
#--------------------------------------------------------------------------------------------------------------------------
# Loading more parameters from the files
L = params['number_of_planes'] # number of planes
material = params['material'] # name of the material or system
face = params['crystal_face'] # crystal face
a0 = params['lattice_parameter'] # lattice parameter of cubic or hexagonal structure
T = params['temperature'] # temperature in K
bc1 = params['BC1_topmost_layer'] # bc value for potential V at the top-most layer
bc2 = params['BC2_in_bulk'] # bc value for V at the bottom-most layer
ef = params['Fermi_level'] # Fermi level in eV as LUL + dE
manifold = params['manifold'] # type of manifold
NWANN = int(params['Wannier_functions_number']) # number of elements in the MLWF basis
HOL = float(params['highest_occupied_level']) # highest occupied level (valence band maximum)
sgeom = params['system_geometry'] # unit-cell symmetry
#--------------------------------------------------------------------------------------------------------------------------
# Creating a logger object for the current program
log_path = identifier + '/' + identifier + '.log'
logger1 = logging.getLogger("BP-eslices")
logger1.setLevel(logging.INFO)
handler = logging.FileHandler(log_path, mode = 'a+')
handler.setFormatter(logging.Formatter("%(levelname)s : %(message)s, line %(lineno)d in %(module)s\n"))
logger1.addHandler(handler)
#--------------------------------------------------------------------------------------------------------------------------
def printlog(text, level = 'i'): # Simple function for general logging
    if level == 'i':
        print(text)
        with open(log_path, mode = 'a+') as OF:
            print(text, file = OF)    
    if level == 'w':
        logger1.warning(text)
        print('WARNING : ' + text + '\n')        
    if level == 'e':
        logger1.error(text)
        raise ValueError(text)  
#--------------------------------------------------------------------------------------------------------------------------
def ChangeBasisRot(xx,yy,fc):
    """ Function to perform the change the basis + rotation algorithm
    according to the confinement direction (fc)."""
    if fc == '100': # the none element is to avoid error when unpacking the returned tuple
        return xx,yy,None 
    if fc == '110':
        Rt = np.array([[1/np.sqrt(2),-1/np.sqrt(2),0.0], # rotation matrix
                       [0.0,0.0,-1.0],
                       [1/np.sqrt(2),1/np.sqrt(2),0.0]])
        b1 = np.array([-0.5,0.5,0.0]) # reciprocal lattice vectors
        b2 = np.array([0.0,0.0,1.0])
        return np.dot(Rt,xx*b1+yy*b2)
    if fc == '111':
        Rt = np.array([[1/np.sqrt(2), -1/np.sqrt(2),0], # rotation matrix
                       [1/np.sqrt(6), 1/np.sqrt(6), -2/np.sqrt(6)],
                       [1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)]])
        b1 = np.array([-1/3,-1/3,2/3]) # reciprocal lattice vectors
        b2 = np.array([2/3,-1/3,-1/3])
        return np.dot(Rt,xx*b1+yy*b2)
    if fc == '001h':
        b1 = np.array([1.0, 1/np.sqrt(3), 0.0]) # reciprocal lattice vectors
        b2 = np.array([0.0, 2/np.sqrt(3), 0.0])
        return xx*b1+yy*b2

#--------------------------------------------------------------------------------------------------------------------------
# MAIN    
#####################################################################################
printlog('\n')
printlog('=========================================================================')
printlog('=========================================================================')
printlog('                                BinPo                                    ')
printlog('=========================================================================')
printlog('=========================================================================')
now = dt.datetime.now()
printlog('\tWELCOME TO THE TIGHT-BINDING POISSON CODE FOR 2DES!')
printlog('\n')
printlog('Author: Emanuel A. Martinez, Universidad Complutense de Madrid, Spain.')
printlog('\n')
printlog('---------------------------------------------------------------------')
printlog('\t\tENERGY SLICES CALCULATION')
printlog('---------------------------------------------------------------------')
printlog('\n')
printlog('BinPo, TB-Poisson Solver for 2DES\n'\
 'Copyright (C) 2021 BinPo Team\n\n'\
 'BP-energy_slices.py is part of BinPo.\n\n'\
 'BinPo is free software: you can redistribute it and/or modify\n'\
 'it under the terms of the GNU General Public License as published by\n'\
 'the Free Software Foundation, either version 3 of the License, or\n'\
 '(at your option) any later version.\n\n'\
 'BinPo is distributed in the hope that it will be useful,\n'\
 'but WITHOUT ANY WARRANTY; without even the implied warranty of\n'\
 'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n'\
 'GNU General Public License for more details.\n\n'\
 'You should have received a copy of the GNU General Public License\n'\
 'along with BinPo. See ~/COPYING file or <https://www.gnu.org/licenses/>.')

if dk_x > 1.0 or dk_y > 1.0:
    printlog('The values for the k-grid shift exceeds the BZ1 periodicity!', label  = 'w')

# Generation of k-grid    
Kmesh = BPM.Kmeshgrid2(Nk, scale = kbox, delta_kx = dk_x, delta_ky = dk_y, a = a0) 

now = dt.datetime.now()
printlog('\n')
printlog('---------------------------------------------------------------------')
printlog('Starting on ' + now.strftime("%d%b%Y") + ' at ' + now.strftime("%H:%M:%S"))
printlog('---------------------------------------------------------------------')
printlog('\n')
printlog('DETAILS :')
printlog('\tIdentifier : ' + identifier)
printlog('\tSurface : ' + material + '(' + face + ')')
printlog('\tNumber of planes : ' + str(L))
printlog('\tK-grid : ' + str(Nk) + ' x ' + str(Nk))
printlog('\t\tK-box factor : '  + str(kbox))
printlog('\t\tK-shift : (' + str(dk_x) + ', ' + str(dk_y) + ')')
printlog('\tNumber of batches : ' + str(Nfr))
printlog('\tTemperature : ' + str(T) + ' K')
printlog('\tFermi level : ' + "{:.5f}".format(ef) + ' eV')
printlog('\tEnergy cut : ' + str(Ecut) + ' eV')
if outfile == 'default':
    printlog('\tOutput file : ' + identifier + '_ES.dat')
else:
    printlog('\tOutput file: ' + outfile + '.dat')
printlog('\n')

start = t.time() # general time reference
#--------------------------------------------------------------------------------------
# files loading and Wannier separation
printlog('Loading files...')
DD  = {} # dictionary to save the r-space Hamiltonian elements separated by planes
for z in os.listdir('./Hr' + material + face):
    Z = np.loadtxt('./Hr' + material + face + '/' + z)
    DD[z.split('.')[0]] = BPM.Wann_Sep(Z, NWANN)
#--------------------------------------------------------------------------------------
Kgrid = np.split(Kmesh,Nfr) # here, the coarse k-grid is split into Nfr batches of size Nk**2/Nfr
                            # in order to reduce the computational cost
V = BPM.PotentialEnergy(L, bc1, BC2 = bc2) # initializing potential instance
V.update_potential(TBP.T[1]) # updating the potential instance with the SC potential values
printlog('Constructing Hv tensor...')
HV = V.to_tensor(tensor_size = int(Nk**2/Nfr), nwann= NWANN) # creating a potential energy tensor
printlog('Done!')
printlog('\n')
printlog("------------------------------------------------------------------------")
    
L0 = [] # auxiliary list to save the data
   
for i in range(Nfr): # main loop in which the energy slice is computed batch by batch 
    printlog('BATCH #'+str(i+1))
    printlog('\n')
    printlog('Hopping calculation and tensors construction...')
    t0 = t.time()
    DDT = {} # dictionary to save the k-space Hamiltonian elements
    for key, val in DD.items(): # creating empty list inside DDT to set the size
                DDT['T_' + key] = []

    for key, val in DD.items(): # filling DDT with the 2D Fourier transformed Hamiltonian elements
                DDT['T_' + key] = BPM.Hopping2D(Kgrid[i], val)
    
    H = BPM.Quasi2DHamiltonian(DDT, L, NWANN) # initializing the Quasi2D Hamiltonian instance
    BPM.Quasi2DHamiltonian.set_parameters(T, ef, len(Kgrid[i]), HOL) # number of k-points is now the contained in each batch
    HK = H.HamiltonianTensor(NWANN) # creating a Hamiltonian tensor
    printlog('Done!')
    printlog('Time spent: ' + "{:.2f}".format(t.time()-t0) + ' seg')
    printlog('\n')      
    printlog('Diagonalizating the Hamiltonian tensor')
    
    # It tries to use scipy routine, which is faster. Otherwise, numpy routine will be used.
    try:
        for j in range(len(Kgrid[i])): 
            # for each k-point in the current batch, the corresponding tensor elements
            # are added and diagonalizing to solve the eigenvalue and eigenvector problem
            # note that the subset of eigenvalues during diagonalization is constrained to the 
            # energy window (wec) under the Ecut. Modify in subset_by_value option the position 
            # of wec respect to the Ecut if you want it above Ecut or with Ecut in the middle of the window.
            HT = HK[j] + HV[j]
            w, v = LA.eigh(HT, lower = True, overwrite_a = True, subset_by_value = (Ecut + ef - wec, Ecut + ef))
            
            x = (Kgrid[i]).T[0][j]
            y = (Kgrid[i]).T[1][j]
            # just if there are eigenvalues within the defined range they will be saved                                 
            if len(w) != 0: # and the eigenvectors will be used to compute the orbital character
                Cxy = BPM.orb_char(v[:,len(w)-1], 'xy', L) # this calculation does not represent a heavy load
                Cyz = BPM.orb_char(v[:,len(w)-1], 'yz', L) # however, at present it makes sense with a t2g manifold
                Czx = BPM.orb_char(v[:,len(w)-1], 'zx', L) 
                #----------------------------------------------------------------
                # Applying the change of basis + rotation algorithm
                xx, yy, _ = ChangeBasisRot(x,y,face)
                #----------------------------------------------------------------
                L0.append([xx,yy,w[len(w)-1],Cxy,Cyz,Czx]) # appending the values of kx, ky, E and projections
    except:
        for j in range(len(Kgrid[i])):
            HT = HK[j] + HV[j]
            w, v = np.linalg.eigh(HT, UPLO = 'L')
    
            x = (Kgrid[i]).T[0][j]
            y = (Kgrid[i]).T[1][j]
            # the eigenvalues are filtered after solving the entire diagonalization in this case                     
            Ev = np.max(np.where((w-ef) < Ecut,(w-ef),-1))
            index = np.where((w-ef) == Ev)
            
            if index[0].size > 0: # this part is a little modified but it does the same as above
                Cxy = 0.0
                Cyz = 0.0
                Czx = 0.0
                for indices in index[0]:
                    Cxy += BPM.orb_char(v[ :,int(indices)], 'xy', L)
                    Cyz += BPM.orb_char(v[ :,int(indices)], 'yz', L)
                    Czx += BPM.orb_char(v[ :,int(indices)], 'zx', L)
                    #----------------------------------------------------------------
                    # Applying the change of basis + rotation algorithm
                    xx, yy, _ = ChangeBasisRot(x,y,face)
                    #----------------------------------------------------------------                    
                L0.append([x,y,Ev,Cxy,Cyz,Czx]) # appending the values of kx, ky, E and projections
                
        printlog('NumPy routine was executed by exception!')        
    
    printlog('\n')
    printlog('Time spent to now: '+"{:.2f}".format(t.time()-start) + '  seg')    
    printlog("------------------------------------------------------------------------")
    
printlog('\n')
printlog('Saving file...')
if outfile == 'default':
    np.savetxt(identifier + '/' + identifier + '_ES.dat', np.array(L0))
else:
    np.savetxt(identifier + '/' + outfile + '.dat', np.array(L0))
    
L0.clear() # cleaning the auxiliary list

end = t.time()
printlog('\n')
printlog('Total time spent: '+"{:.2f}".format(end-start) + '  seg') 
printlog('\n')
printlog('============================================================================')
now2 = dt.datetime.now()
printlog('Finishing on ' + now2.strftime("%d%b%Y") + ' at ' + now2.strftime("%H:%M:%S"))
printlog('CALCULATION DONE!')
printlog('============================================================================')
printlog('\n')

