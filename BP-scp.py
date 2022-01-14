#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 BinPo, TB-Poisson Solver for 2DES
 Copyright (C) 2021 BinPo Team
 
 BP-scp.py is part of BinPo.
 
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
     SELF-CONSISTENT POTENTIAL ENERGY (V) CALCULATION
     BP-scp.py component perform the TB-Poisson scheme in order to obtain 
     the self-consistent potential of the system.
"""

__author__ = "Emanuel A. Martinez"
__email__ = "emanuelm@ucm.es"
__copyright__ = "Copyright (C) 2021 BinPo Team"
__version__ = 1.0
__date__ = "January 14, 2022"

import os
import numpy as np
import matplotlib.pyplot as plt
import BPmodule as BPM
import BPdatabase as BPD
import logging
import time as t
import datetime as dt
import argparse
import yaml

# Loading the configuration file to set parameters not defined by terminal
#--------------------------------------------------------------------------------------------------------------------------
with open('./config_files/scp.yaml','r') as f:
    try:
        data = yaml.load(f, Loader=yaml.FullLoader)
    except yaml.YAMLError as exc:
        print ("Error in configuration file: ", exc)
#--------------------------------------------------------------------------------------------------------------------------
# Messages to print by terminal when checking the description with -h or --help command
description = "SELF-CONSISTENT POTENTIAL ENERGY (V) CALCULATION"
epilog = "By default, arguments not specified are taken from './config_files/scp.yaml'."
 
msj_id = "Identifier for the calculation. String. It is unique and can be recalled later in post-processing."
msj_mat = "Name of the material. String. It could be 'STO' or 'KTO'."
msj_fc = "Confinement direction 'hkl'. String. Allowed values: '100' ('010','001'), '110' ('011','101'), '111'."
msj_T = "Temperature in K. Float. Allowed range: [0.001, 315.0]"
msj_L = "Number of planes stacked along normal direction. Integer. Allowed range: [10, 200]."
msj_Nk = "Square root of the total number of points in 2D k-grid. Integer. It must be > 10."
msj_BC1 = "Dirichlet boundary condition for V in eV at the top-most layer (i.e. V[0]). Float. Allowed range: [-1, -0.001]."
msj_BC2 = "Dirichlet boundary condition for V in eV in the bottom-most layer (i.e. V[L-1]). Float. Allowed range: [bc1, 3*|bc1|]."
msj_dE = "Energy shift (dE) from LUL in eV. Used to define the Fermi level as ef = LUL + dE. Float. Allowed range: [0.0, 0.02]."
msj_pm = "Model used for relative permittivity as function of electric field E. String. Available models are: 'Cop', 'Nev1', 'Nev2', 'Mat' and 'cte_N'. "\
"Alternatively, an expression as function of E in python syntax could be introduced. See the comments in scp.yaml for information."
msj_mf = "Mixing factor for the over-relaxation mixing method. Float. Tipically 0.06-0.4."
#--------------------------------------------------------------------------------------------------------------------------
# Passing arguments by terminal and setting their default values
parser = argparse.ArgumentParser(description = description, epilog = epilog)
parser.add_argument('-id', default = data['SCP_CALCULATION']['identifier'], type = str, help = msj_id)
parser.add_argument('-mt', default = data['SCP_CALCULATION']['material'], type = str, help = msj_mat)
parser.add_argument('-cfd', default = data['SCP_CALCULATION']['crystal_face'], type = str, help = msj_fc)
parser.add_argument('-te', default = data['SCP_CALCULATION']['temperature'], type = float, help = msj_T)
parser.add_argument('-tl', default = data['SCP_CALCULATION']['number_of_planes'], type = int, help = msj_L)
parser.add_argument('-nk', default = data['SCP_CALCULATION']['sqrt_kgrid_numbers'], type = int, help = msj_Nk)
parser.add_argument('-bc1', default = data['SCP_CALCULATION']['BC1_topmost_layer'], type = float, help = msj_BC1)
parser.add_argument('-bc2', default = data['SCP_CALCULATION']['BC2_in_bulk'], type = float, help = msj_BC2)
parser.add_argument('-de', default = data['SCP_CALCULATION']['shift_from_LUL'], type = float, help = msj_dE)
parser.add_argument('-pm', default = data['SCP_CALCULATION']['permittivity_model'], type = str, help = msj_pm)
parser.add_argument('-mf', default = data['SCP_CALCULATION']['mixing_factor'], type = float, help = msj_mf)
args = parser.parse_args()
#--------------------------------------------------------------------------------------------------------------------------
# Setting at once what depends on the identifier
identifier = args.id # Identifier for the calculation
# Checking if a directory with the same identifier name already exists! If not, it creates one
if os.path.isdir(identifier) == True: 
    raise ValueError("identifier '%s' already exists! Please, remove the folder or use another identifier." % identifier)
os.system('mkdir ' + identifier)
#--------------------------------------------------------------------------------------------------------------------------
# Creating a logger object for the current program
log_path = identifier + '/' + identifier + '.log'
logger1 = logging.getLogger("BP-scp")
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
# Copying the rest of parameters updatable by terminal
#--------------------------------------------------------------------------------------------------------------------------
material = (args.mt).upper() # uppercase the material name
face_in = args.cfd # crystal face
T = args.te # temperature in K
L = args.tl # number of planes
Nk = args.nk # (total number of kpoints = Nk*Nk)
bc1 = args.bc1 # boundary condition at z = 0 (i.e. at the top-most layer)
bc2 = args.bc2 # boundary condition at z = L-1 (i. e. bulk)
dE = args.de #shift from Fermi level
permitt_model = args.pm #permittivity model
Fmix =  args.mf# mixing factor
#---------------------------------------------------------------------
# These parameters are loaded from conf_preproc.yaml file!!
filename = BPD.materials[material]['Wannier_file']
a0 = BPD.materials[material]['lattice_parameter'] # lattice parameter of cubic structure
LUL = BPD.materials[material]['lowest_unoccupied_level'] # lowest unoccupied level
#---------------------------------------------------------------------
# These parameters are loaded from scp_preproc.yaml file!!
visualization = data['SCP_CALCULATION']['potential_live_visualization'] # whether or not live visualization of V(z)
err_visualization = data['SCP_CALCULATION']['error_live_visualization'] # whether or not live visualization of error vs iterations

conv_thr = data['SCP_CALCULATION']['ADVANCED_PARAMETERS']['conv_threshold'] # convergence value when comparing potentials
MAX_ITER = data['SCP_CALCULATION']['ADVANCED_PARAMETERS']['max_iterations'] # total allowed iterations for TB-Poisson scheme
Vinitial = data['SCP_CALCULATION']['ADVANCED_PARAMETERS']['V_initial'] # Initial shape of the potential energy

method = data['SCP_CALCULATION']['ADVANCED_PARAMETERS']['Total_Hk_method'] # method to create the Hamiltonian tensor
    
Neumann_condition = data['SCP_CALCULATION']['Neumann_at_bulk'] # use or not Neumann b. c. at bottom-most layer
k_shift = data['SCP_CALCULATION']['k_shift'] # Displacements of the k-grid

#####################################################################################
# In the case that a fixed background charge is used
use_fixQ = data['SCP_CALCULATION']['ADVANCED_PARAMETERS']['cons_charge_background']
Qdef = data['SCP_CALCULATION']['ADVANCED_PARAMETERS']['charge_per_site']
Qextent = data['SCP_CALCULATION']['ADVANCED_PARAMETERS']['charge_extension']
#####################################################################################
#--------------------------------------------------------------------------------------------------------------------------
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
# Checking one by one if the values for the input parameters are correct
if material != 'STO' and material != 'KTO':
    printlog("%s is and invalid material name!", level = 'e')
if face_in not in ['100','010','001','110','101','011','111']:
    printlog('%s are invalid Miller indices' % face_in, level = 'e')
if T < 1.0e-3 or T > 315.0 :
    printlog('temperature is out of allowed range!', level = 'e')
if a0 <= 0:
    printlog('lattice parameter must be positive!', level = 'e')  
if L <= 10 or L > 200 :
    printlog('number_of_planes must range from 10 to 200!', level = 'e')
if Nk < 10:
    printlog('sqrt_kgrid_numbers must be greater than 10!', level = 'e')
if bc1 >= -1e-3 or bc1 <= -1.0:
    printlog('BC1_topmost_layer must range from -1 to -0.001 eV', level = 'e')
if bc2 < (bc1 + 0.01*abs(bc1)) or bc2 > 3*abs(bc1):
    printlog('BC2_in_bulk must range from BC1 to 3*|BC1|!', level = 'e')
if dE < 0.0 or dE > 0.02:
    printlog('shift_from_LUL must be between 0 and 0.02 eV', level = 'e')

eps_list = ['Cop', 'Ang', None]
if permitt_model not in eps_list and 'cte' not in permitt_model: # Check if permittivity is in the list of available ones
    E = 100 # test value for checking if permittivity expression is correct.
    try:
        eval(permitt_model)
    except:
        printlog('%s is invalid model or expression! Please, make sure of either selecting a permittivity'\
                         ' from the available ones or correctly using Python syntax!' % permitt_model, level = 'e')

if method != 'vectorized' and method != 'iterable':
    printlog("%s is an invalid method. Available methods are 'vectorized' and 'iterable'." % method, level = 'e')

if use_fixQ == True:
    if Qextent < 0 or Qextent > L:
        printlog('Fixed density background out of valid range!', level = 'e')
    if Qdef >= 0.05:
        printlog('The fixed charge per site seems to be very high!', level = 'w')
# Creating specific quantities directly from parameters
ef = LUL + dE # Fermi level in eV
Kmesh = BPM.Kmeshgrid(Nk, delta_kx = k_shift[0], delta_ky = k_shift[1]) # Generation of k-grid
crystal = BPM.CrystalFeatures(face_in, a0, material) # Initialization of crystal properties
face = crystal.face # It turns the indices into an unique ones per face ['100', '110', '111']

printlog('---------------------------------------------------------------------')
printlog('\tCALCULATION OF SELF-CONSISTENT POTENTIAL ENERGY')
printlog('---------------------------------------------------------------------')
printlog('\n')
printlog('BinPo, TB-Poisson Solver for 2DES\n'\
 'Copyright (C) 2021 BinPo Team\n\n'\
 'BP-scp.py is part of BinPo.\n\n'\
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
printlog('\n')
printlog('---------------------------------------------------------------------')
printlog('\tStarting on ' + now.strftime("%d%b%Y") + ' at ' + now.strftime("%H:%M:%S"))
printlog('---------------------------------------------------------------------')
printlog('\n')
printlog('DETAILS:')
printlog('\tIdentifier : ' + identifier)
printlog('\tSurface : ' + material + '(' + face + ')')
printlog('\tNumber of planes : ' + str(L))
printlog('\tK-grid : ' + str(Nk) + ' x ' + str(Nk))
printlog('\t\tK-shift : (' + str(k_shift[0]) + ' ,' + str(k_shift[1]) + ')')
printlog('\tBoundary conditions : ')
printlog('\t\tV[0] = ' + str(bc1) + ' eV')
printlog('\t\tV[L-1] = ' + str(bc2) + ' eV')
printlog('\tNeumann condition at V[L-1] : ' + str(Neumann_condition))
printlog('\tPermittivity model : ' + permitt_model)
printlog('\tTemperature : ' + str(T) + ' K')
printlog('\tFermi level : ' + "{:.5f}".format(ef) + ' eV')
printlog('\tTotal Hk method : ' + method)
printlog('\tMixing factor : ' + str(Fmix))
printlog('\tConvergence threshold : ' + str(conv_thr))
printlog('\tUsing charge background : ' + str(use_fixQ))
if use_fixQ == True:
    printlog('\t\tCharge per site : ' + str(Qdef))
    printlog('\t\tCharge extent : ' + str(Qextent))
printlog('\tInitial V shape : ' + Vinitial)
printlog('\n')
start = t.time() # general time reference
if Fmix < 0.05:
    printlog("A mixing factor of %s could be too small. The total number of iterations could considerably increase" % Fmix, level = 'w')
if Fmix > 0.4:
    printlog("A mixing factor of %s could be too big. It could result in a bouncing potential without reaching convergence"\
             " or requiring a large number of iterations" % Fmix, level = 'w')
#--------------------------------------------------------------------------------------
# files loading and Wannier separation
printlog('Loading files...')
Z_7 = np.loadtxt('./Hr' + material + face + '/Z_7.dat')
Z_6 = np.loadtxt('./Hr' + material + face + '/Z_6.dat')
Z_5 = np.loadtxt('./Hr' + material + face + '/Z_5.dat')
Z_4 = np.loadtxt('./Hr' + material + face + '/Z_4.dat')
Z_3 = np.loadtxt('./Hr' + material + face + '/Z_3.dat')
Z_2 = np.loadtxt('./Hr' + material + face + '/Z_2.dat')
Z_1 = np.loadtxt('./Hr' + material + face + '/Z_1.dat')
Z0 = np.loadtxt('./Hr' + material + face + '/Z0.dat')

D_7 = BPM.Wann_Sep(Z_7) # Separation of the <Pw|H|Pw'> for the plane P 
D_6 = BPM.Wann_Sep(Z_6) # and the Wannier functions w, w'
D_5 = BPM.Wann_Sep(Z_5)
D_4 = BPM.Wann_Sep(Z_4)
D_3 = BPM.Wann_Sep(Z_3)
D_2 = BPM.Wann_Sep(Z_2)
D_1 = BPM.Wann_Sep(Z_1)
D0 = BPM.Wann_Sep(Z0)
#--------------------------------------------------------------------------------------
# 2D Fourier transform of the <Pw|H|Pw'> elements
printlog("Transforming <0w|H|Rw'> to k-space...")
t0 = t.time()
T_7 = BPM.Hopping2D(Kmesh, D_7)
T_6 = BPM.Hopping2D(Kmesh, D_6)
T_5 = BPM.Hopping2D(Kmesh, D_5)
T_4 = BPM.Hopping2D(Kmesh, D_4)
T_3 = BPM.Hopping2D(Kmesh, D_3)
T_2 = BPM.Hopping2D(Kmesh, D_2)
T_1 = BPM.Hopping2D(Kmesh, D_1)
T0 = BPM.Hopping2D(Kmesh, D0)
#--------------------------------------------------------------------------------------
# Initializing the Hamiltonian
BPM.Quasi2DHamiltonian.set_parameters(T, ef, Nk*Nk)
H = BPM.Quasi2DHamiltonian(T_7, T_6, T_5, T_4, T_3, T_2, T_1, T0, L)
printlog('Done!')
printlog('Time spent: ' + "{:.2f}".format(t.time()-t0) + ' seg')
printlog('\n')
#--------------------------------------------------------------------------------------
t1 = t.time()
# If the method to generate the Hamiltonian H is 'vectorized' the whole H is created at once
if method == 'vectorized':
    printlog('Constructing the slab Hamiltonian...')
    Hk = H.HamiltonianTensor()
    printlog('Done!')
    printlog('Time spent: ' + "{:.2f}".format(t.time()-t1) + ' seg')
    printlog('\n')
#--------------------------------------------------------------------------------------
# Here we create the potential instance and set the crystal properties
printlog('Setting crystal properties and initializing potential energy...')
V = BPM.PotentialEnergy(L, bc1, bc2, permitt_model, V_init = Vinitial)
d_hkl = crystal.interplanar_distance() 
A_hkl = crystal.face_area()
printlog('Done!')
printlog('\n')
#--------------------------------------------------------------------------------------
# whether or not to use the live visualizations
if visualization == True:
    plt.style.use('ggplot')
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(1,1,1)
    plt.subplots_adjust(left = 0.14, bottom = 0.14)
    
if err_visualization == True:
    plt.style.use('ggplot')
    fig2 = plt.figure(figsize = (3,3))
    ax2 = fig2.add_subplot(1,1,1)
    plt.subplots_adjust(left = 0.2)    
#-------------------------------------------------------------------------------------- 
# Starting the Tight-Binding Poisson iterations
err = 1 # We choose an arbitrary initial value for total error
Lerr = [] # Into this list the errors will be append if err_visualization == True.
while err > conv_thr:
    printlog('============================================================================')
    printlog('\n')
    printlog('-------------------')
    printlog('ITERATION #'+str(V.counter))
    printlog('-------------------')
    printlog('\n')
    printlog("------------------------------------Vi------------------------------------")
    printlog(V)
    printlog('\n')
    printlog('-------------------')
    printlog('Tight-binding step:')
    printlog('-------------------')
    printlog('\n')
    t2 = t.time()
    printlog('Starting diagonalization and charge density calculation...')
    # Charge solver is computed according to the method selected
    if method == 'vectorized':
        rho = BPM.Quasi2DHamiltonian.ChargeSolver(Hk, V)
    if method == 'iterable':
        rho = BPM.Quasi2DHamiltonian.ChargeSolver2(H, V)      
    printlog('Done!')
    printlog('Time spent: ' + "{:.2f}".format(t.time()-t2) + ' seg')
    printlog('\n')
    printlog('-------------')
    printlog('Poisson step:')
    printlog('-------------') 
    printlog('\n')
    t3 = t.time()
    printlog('Solving Poisson equation...')
    Vout, n_iter = V.Poisson_solver(rho ,d_hkl ,A_hkl, neumann_bc2 = Neumann_condition, background_Q = use_fixQ, site_qdef = Qdef, ext_qdef = Qextent)
    printlog('Iteration converged after ' + str(n_iter) + ' steps.')
    printlog('Time spent: ' + "{:.2f}".format(t.time()-t3) + ' seg')
    printlog('\n')
    err = V.Error(Vout)
    Vin = np.copy(V.values)
    printlog("------------------------------------Vout------------------------------------")
    printlog('V = ')
    printlog(Vout)
    printlog('\n')
    printlog("TOTAL ERROR = "+ str(err))
    delta_t = t.time()-start
    printlog("Time spent to now: " + "{:.2f}".format(delta_t) + '  seg')
    printlog('\n')
    printlog('\n')
    
    if visualization == True:
        BPM.SnapPlotter(fig1, ax1, V, Vin, Vout, delta_t, err, 0.5)
        
    if err_visualization == True:
        Lerr.append(err)
        BPM.SnapError(fig2, ax2, Lerr, conv_thr, 1.0)
    
    if err < conv_thr: # If SC is achieved the following quantities are saved 
        V.update_potential(Vout)
        Lout = []
        Lout.append(np.arange(V.L)) # planes
        Lout.append(V.values) # SC potential energy
        Lout.append(rho) # electron density in number of electrons/plane
        Lout.append(rho/(A_hkl*d_hkl*1e6)) # electron density in number of electrons/cm**2
        Lout.append(V.electric_field(d_hkl)) # electric field in V/m
        Lout.append(V.permittivity_model(V.electric_field(d_hkl))) # relative permittivity
        printlog('\n')
        printlog('Saving files...\n')
        np.savetxt(identifier + '/' + identifier + '_SCP.dat', np.array(Lout).T)
        # After saving the SC-solution the parameters of this problem are also saved!
        # Creating a dictionary to make a .YAML file for post-processing access to the parameters
        dict_file = {'identifier': identifier,
                     'material' : material,
                     'lattice_parameter': a0,
                     'crystal_face' : face,
                     'number_of_planes' : L,
                     'sqrt_kgrid_numbers' : Nk,
                     'k_shift' : k_shift,
                     'temperature' : T,
                     'Fermi_level' : ef,
                     'BC1_topmost_layer' : bc1,
                     'BC2_in_bulk' : bc2,
                     }
        with open(identifier + '/' + identifier + '.yaml', 'w') as file:
            documents = yaml.dump(dict_file, file)
        break
    
    #==============================================================
    # Mixing the potentials with the over-relaxation method
    Vnew = Vin + Fmix*(Vout - Vin) 
    V.update_potential(Vnew) # Update the PotentialEnergy intance  
    #==============================================================
    
    if V.counter == MAX_ITER:
        printlog('\n')    
        end = t.time()
        printlog('============================================================================')
        printlog('Convergence not achieved after ' + str(V.counter) + ' iterations.\nTotal time spent: ' + "{:.2f}".format(end-start) + '  seg')
        printlog('\n')
        now2 = dt.datetime.now()
        printlog('Finishing on ' + now2.strftime("%d%b%Y") + ' at ' + now2.strftime("%H:%M:%S"))
        printlog('============================================================================')
        printlog('Runtime Error : Convergence not achieved after ' + str(int(MAX_ITER)) + ' iterations.')
        raise RuntimeError('Convergence not achieved after ' + str(int(MAX_ITER)) + ' iterations.')

   
printlog('\n')    
end = t.time()
printlog('============================================================================')
printlog('Convergence achieved after ' + str(V.counter) + ' iterations.\nTotal time spent: ' + "{:.2f}".format(end-start) + '  seg')
printlog('\n')
now2 = dt.datetime.now()
printlog('Finishing on ' + now2.strftime("%d%b%Y") + ' at ' + now2.strftime("%H:%M:%S"))
printlog('CALCULATION DONE!')
printlog('============================================================================')
printlog('\n')

