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
     bandstructure.
"""
import numpy as np
import BPmodule as BPM
import time as t
import datetime as dt
import argparse
import yaml
import scipy.linalg as LA
import logging

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
# Copying the rest of the parameters updatable by terminal
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
material = params['material']
face = params['crystal_face'] # crystal face
a0 = params['lattice_parameter'] # lattice parameter of cubic structure 
T = params['temperature'] # temperature in K
bc1 = params['BC1_topmost_layer']
bc2 = params['BC2_in_bulk']
ef = params['Fermi_level'] #Fermi level in eV as LUL + dE
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
Kmesh = BPM.Kmeshgrid(Nk, scale = kbox, delta_kx = dk_x, delta_ky = dk_y) # Generation of kgrid
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

printlog('Loading files...')
Z_7 = np.loadtxt('./Hr' + material + face + '/Z_7.dat')
Z_6 = np.loadtxt('./Hr' + material + face + '/Z_6.dat')
Z_5 = np.loadtxt('./Hr' + material + face + '/Z_5.dat')
Z_4 = np.loadtxt('./Hr' + material + face + '/Z_4.dat')
Z_3 = np.loadtxt('./Hr' + material + face + '/Z_3.dat')
Z_2 = np.loadtxt('./Hr' + material + face + '/Z_2.dat')
Z_1 = np.loadtxt('./Hr' + material + face + '/Z_1.dat')
Z0 = np.loadtxt('./Hr' + material + face + '/Z0.dat')
    
D_7 = BPM.Wann_Sep(Z_7)    
D_6 = BPM.Wann_Sep(Z_6)
D_5 = BPM.Wann_Sep(Z_5)
D_4 = BPM.Wann_Sep(Z_4)
D_3 = BPM.Wann_Sep(Z_3)
D_2 = BPM.Wann_Sep(Z_2)
D_1 = BPM.Wann_Sep(Z_1)
D0 = BPM.Wann_Sep(Z0)
 
Kgrid = np.split(Kmesh,Nfr)
    
V = BPM.PotentialEnergy(L, bc1, BC2 = bc2)
V.update_potential(TBP.T[1])
printlog('Constructing Hv tensor...')
HV = V.to_tensor(int(Nk**2/Nfr))
printlog('Done!')
    
L0 = []# to save the data
    
for i in range(Nfr):
    printlog('BATCH #'+str(i+1))
    printlog('\n')
    printlog('Hopping calculation and tensors construction...')
    t0 = t.time()
    
    T_7 = BPM.Hopping2D(Kgrid[i], D_7)    
    T_6 = BPM.Hopping2D(Kgrid[i], D_6)
    T_5 = BPM.Hopping2D(Kgrid[i], D_5)
    T_4 = BPM.Hopping2D(Kgrid[i], D_4)
    T_3 = BPM.Hopping2D(Kgrid[i], D_3)
    T_2 = BPM.Hopping2D(Kgrid[i], D_2)
    T_1 = BPM.Hopping2D(Kgrid[i], D_1)
    T0 = BPM.Hopping2D(Kgrid[i], D0)
        
    printlog("Transforming <0w|H|Rw'> to k-space...")    
    BPM.Quasi2DHamiltonian.set_parameters(T,ef,len(Kgrid[i]))
    H = BPM.Quasi2DHamiltonian(T_7, T_6, T_5, T_4, T_3, T_2, T_1, T0, L)
    
    HK = H.HamiltonianTensor()
    printlog('Done!')
    printlog('Time spent: '+"{:.2f}".format(t.time()-t0) + '  seg')
    printlog('\n')
        
    printlog('Diagonalizating the Hamiltonian tensor')
    t0 = t.time()
    
    try:
        for j in range(len(Kgrid[i])):
            HT = HK[j] + HV[j]
            w, v = LA.eigh(HT, lower = True, overwrite_a = True, subset_by_value = (Ecut + ef - wec, Ecut + ef))
            
            x = (Kgrid[i]).T[0][j]
            y = (Kgrid[i]).T[1][j]
                                             
            if len(w) != 0:
                Cxy = BPM.orb_char(v[:,len(w)-1], 'xy', L)
                Cyz = BPM.orb_char(v[:,len(w)-1], 'yz', L)
                Czx = BPM.orb_char(v[:,len(w)-1], 'zx', L)
                
                L0.append([x,y,w[len(w)-1],Cxy,Cyz,Czx])
    except:
        for j in range(len(Kgrid[i])):
            HT = HK[j] + HV[j]
            w, v = np.linalg.eigh(HT, UPLO = 'L')
    
            x = (Kgrid[i]).T[0][j]
            y = (Kgrid[i]).T[1][j]
                                 
            Ev = np.max(np.where((w-ef) < Ecut,(w-ef),-1))
            index = np.where((w-ef) == Ev)
            
            if index[0].size > 0:
                Cxy = 0.0
                Cyz = 0.0
                Czx = 0.0
                for indices in index[0]:
                    Cxy += BPM.orb_char(v[ :,int(indices)], 'xy', L)
                    Cyz += BPM.orb_char(v[ :,int(indices)], 'yz', L)
                    Czx += BPM.orb_char(v[ :,int(indices)], 'zx', L)
                L0.append([x,y,Ev,Cxy,Cyz,Czx])
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

