# -*- coding: utf-8 -*-
"""
 BinPo, TB-Poisson Solver for 2DES
 Copyright (C) 2021 BinPo Team
 
 BP-orb_density.py is part of BinPo.
 
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
     ORBITAL DECOMPOSITION OF ELECTRON DENSITY
     BP-orb_density.py component allows for decomposing and plotting the 
     electron density according to the orbital character. At present, this
     component is only available for cubic systems with t2g manifold, such
     as STO, KTO and CTO.
"""

__author__ = "Emanuel A. Martinez"
__email__ = "emanuelm@ucm.es"
__copyright__ = "Copyright (C) 2021 BinPo Team"
__version__ = 1.1
__date__ = "August 9, 2022"

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as LA
import BPmodule as BPM
import time as t
import datetime as dt
import argparse
import yaml
import logging
import os

# Loading the configuration files to set parameters not defined by terminal
#--------------------------------------------------------------------------------------------------------------------------
with open('./config_files/orb_density.yaml', 'r') as f:
    try:
        data = yaml.load(f, Loader = yaml.FullLoader)
    except yaml.YAMLError as exc:
        print ("Error in configuration file: ", exc)
#--------------------------------------------------------------------------------------------------------------------------
# Messages to print by terminal when checking the description with -h or --help command
msj_id = 'Identifier to compute the orbital decomposition of the electron density. String.'
msj_Nk = "Square root of the total number of points in 2D k-grid. Integer. It must be > 10."
msj_xy = 'Limits for the plot. Float (list). Introduce it as x_min x_max y_min y_max. If not, the limits will be authomatically set.'
msj_aa = 'Opacity of the curves shadowing. Float. It must be between 0.0 (transparent) and 1.0 (solid).'
description = "ORBITAL DECOMPOSITION OF ELECTRON DENSITY"
epilog = "By default, arguments not specified are taken from './config_files/orb_density.yaml'."
#--------------------------------------------------------------------------------------------------------------------------
# Passing arguments by terminal and setting their default values
parser = argparse.ArgumentParser(description = description, epilog = epilog)
parser.add_argument('-id', type = str, default = data['DENSITY_ANALYSIS']['identifier'], help = msj_id)
parser.add_argument('-nk', type = int, default = data['DENSITY_ANALYSIS']['sqrt_kgrid_numbers'], help = msj_Nk)
parser.add_argument('-aa', type = float, default = data['DENSITY_ANALYSIS']['PLOT_ADJUST']['alpha'], help = msj_aa)
parser.add_argument('-xy', nargs = 4, type = float, help = 'Limits for the plot.')
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
# Copying the rest of the parameters updatable by terminal, except xy
ALPHA = args.aa
Nk = args.nk
#--------------------------------------------------------------------------------------------------------------------------
# Loading more parameters from the files
a0 = params['lattice_parameter'] # lattice parameter of cubic or hexagonal structure 
T = params['temperature'] # temperature in K
bc1 = params['BC1_topmost_layer'] # bc value for potential V at the top-most layer
bc2 = params['BC2_in_bulk'] # bc value for V at the bottom-most layer
ef = params['Fermi_level'] # Fermi level in eV as LUL + dE
L = params['number_of_planes'] # number of planes
material = params['material'] # name of the material or system
face = params['crystal_face'] # confinement direction
dk_x, dk_y = params['k_shift'] # shift in k-grid
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
def PartialChargePlotter(rho0, rho1, rho2, rho3, data3):
    # unfolding data0 dictionary to load several details
    # for the matplotlib plot
    plotstyle = data3['PLOT_ADJUST']['plotstyle']
    fig_size = data3['PLOT_ADJUST']['fig_size']
    lw = data3['PLOT_ADJUST']['linewidth']
    c0, c1, c2, c3 = data3['PLOT_ADJUST']['color_seq'].split(',')
    ms = data3['PLOT_ADJUST']['markersize']
    l3, b3, r3, t3 = data3['PLOT_ADJUST']['axis_adjust']
    ptitle = data3['PLOT_ADJUST']['title']
    ptitle_size = data3['PLOT_ADJUST']['title_size']
    use_fill = data3['PLOT_ADJUST']['fill_curves']
    
    xlabel = data3['LABELS']['xlabel']
    xfontsize = data3['LABELS']['xfontsize']
    ylabel = data3['LABELS']['ylabel']
    yfontsize = data3['LABELS']['yfontsize']
    tick_size = data3['LABELS']['ticksize']
    use_legends = data3['LABELS']['legends']
    legend_size = data3['LABELS']['legend_size']
    
    save_data = data3['SAVING']['save_data']
    save_plot = data3['SAVING']['save_plot']
    pformat = data3['SAVING']['format']
    resol = data3['SAVING']['dpi']
    
    plt.style.use(plotstyle)
    z = np.arange(L)
    d = np.zeros_like(z)
    fig, ax = plt.subplots(figsize = fig_size)
    ax.plot(z, rho0, lw = lw, marker = 'o', ms = ms, label = 'total', color = c0, zorder = 10)
    ax.plot(z, rho1, lw = lw, marker = 's', ms = ms, label = 'xy', color = c1, zorder = 20)
    ax.plot(z, rho2, lw = lw, marker = 'd', ms = ms, label = 'yz', color = c2, zorder = 30)
    ax.plot(z, rho3, lw = lw, marker = '<', ms = ms, label = 'zx', color = c3, zorder = 40)
    if use_fill == True:
        try:
            ax.fill_between(list(z), list(rho0), where = rho0 >= d, interpolate = True, color = c0, alpha = ALPHA, zorder = 10)
            ax.fill_between(list(z), list(rho1), where = rho1 >= d, interpolate = True, color = c1, alpha = ALPHA, zorder = 20)
            ax.fill_between(list(z), list(rho2), where = rho2 >= d, interpolate = True, color = c2, alpha = ALPHA, zorder = 30)
            ax.fill_between(list(z), list(rho3), where = rho3 >= d, interpolate = True, color = c3, alpha = ALPHA, zorder = 40)
        except:
            printlog('The fill under the curves could not be applied!', level = 'w')
    #ax.plot(z, rho1 + rho2 + rho3, lw = lw, ls = 'dotted', label = 'cumulative') # check if the cumulative == total
    ax.set_xlabel(xlabel, size = xfontsize)
    ax.set_ylabel(ylabel, size = yfontsize)
    ax.tick_params(labelsize = tick_size)
    if len(ptitle) > 0:
        ax.set_title(ptitle, fontsize = ptitle_size)
    plt.subplots_adjust(left = l3, bottom = b3, right = r3, top = t3)
    if args.xy:
        ax.set_xlim(args.xy[0], args.xy[1])
        ax.set_ylim(args.xy[2], args.xy[3])
    plt.subplots_adjust(left = 0.15)
    if use_legends == True:
        plt.legend(loc = 'best', facecolor = 'w', framealpha = 0.6, fontsize = legend_size)
    Lout = [z, rho0, rho1, rho2, rho3]
    if save_data == True:
        printlog('Saving data...')
        np.savetxt(identifier + '/' + identifier + '_orb-dens.dat', np.array(Lout).T)    
    if save_plot == True:
    	printlog('Saving plot...')
    	plt.savefig(identifier + '/' + identifier + '_orb-dens' + pformat, dpi = resol)    	
    return None
#--------------------------------------------------------------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------------------------
# NOTE : At moment, this component is only available for cubic systems
# with t2g manifold, such as STO, KTO and CTO.    
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
printlog('\t   ORBITAL DECOMPOSITION OF ELECTRON DENSITY')
printlog('---------------------------------------------------------------------')
printlog('\n')
printlog('BinPo, TB-Poisson Solver for 2DES\n'\
 'Copyright (C) 2021 BinPo Team\n\n'\
 'BP-orb_density.py is part of BinPo.\n\n'\
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
    
Kmesh = BPM.Kmeshgrid(Nk, delta_kx = dk_x, delta_ky = dk_y) # Generation of k-grid

Nk_min = 10 # lower limit for Nk
if Nk < Nk_min:
    printlog('sqrt_kgrid_numbers must be greater than ' + str(Nk_min), level = 'e')

if manifold != 't2g':
    printlog('\n')
    printlog('Orbital decomposition of the electron density onto an arbitrary set of MLWFs is not implemented yet.', 'e')
    printlog('\n')

now = dt.datetime.now()
printlog('\n')
printlog('---------------------------------------------------------------------')
printlog('Starting on ' + now.strftime("%d%b%Y") + ' at ' + now.strftime("%H:%M:%S"))
printlog('---------------------------------------------------------------------')
printlog('\n')
printlog('DETAILS:')
printlog('\tIdentifier: ' + identifier)
printlog('\tSurface : ' + material + '(' + face + ')')
printlog('\tNumber of planes: ' + str(L))
printlog('\tK-grid: ' + str(Nk) + ' x ' + str(Nk))
printlog('\tTemperature: ' + str(T) + ' K')
printlog('\tFermi level: ' + "{:.5f}".format(ef) + ' eV')
printlog('\n')

start = t.time() # general time reference
#############################################################################################################
# files loading and Wannier separation
printlog('Loading files...')
DD  = {} # dictionary to save the r-space Hamiltonian elements separated by planes
for z in os.listdir('./Hr' + material + face):
    Z = np.loadtxt('./Hr' + material + face + '/' + z)
    DD[z.split('.')[0]] = BPM.Wann_Sep(Z, NWANN)

printlog("Transforming <0w|H|Rw'> to k-space...") # 2D Fourier transform of the <0 w|H|Rxy+Rz w'> elements
printlog("It could take a while...")
t0 = t.time()
DDT = {} # dictionary to save the k-space Hamiltonian elements
for key, val in DD.items(): # creating empty list inside DDT to set the size
            DDT['T_' + key] = []

for key, val in DD.items(): # filling DDT with the 2D Fourier transformed Hamiltonian elements
            DDT['T_' + key] = BPM.Hopping2D(Kmesh, val)

H = BPM.Quasi2DHamiltonian(DDT, L, NWANN) # initializing the Quasi2D Hamiltonian instance
BPM.Quasi2DHamiltonian.set_parameters(T, ef, Nk*Nk, HOL)
printlog('Done!')
printlog('Time spent: ' + "{:.2f}".format(t.time()-t0) + ' seg')
printlog('\n')
#############################################################################################################

t1 = t.time()
printlog('Constructing the slab Hamiltonian...')
Hk = H.HamiltonianTensor(NWANN) # creating Hamiltonian tensor
printlog('Done!')
printlog('Time spent: ' + "{:.2f}".format(t.time()-t1) + ' seg')

printlog('Setting crystal properties and initializing potential energy...')
V = BPM.PotentialEnergy(L, bc1, bc2) # initializing potential instance
V.update_potential(TBP.T[1]) # updating the potential instance with the SC potential values
crystal = BPM.CrystalFeatures(face, a0, material) # initializing the crystal instance
d_hkl = crystal.interplanar_distance() # interplanar distance
A_hkl = crystal.face_area() # area of 2D real lattice
printlog('Done!')

t2 = t.time()
printlog('Starting diagonalization and charge density calculation...')
printlog('It could take a while...')

Hv = V.to_tensor() # creating potential energy tensor
Lrho = [] # auxiliary list to save the electron density
Lrho_xy, Lrho_yz, Lrho_zx = [], [], [] # auxiliary list to save the electron density
                                       # according to the orbital character 
#--------to select the orbital indices
zx_u = np.arange(0,6*L,6)
zx_d = np.arange(1,6*L,6)
zx = np.sort(np.concatenate([zx_u,zx_d]))
yz_u = np.arange(2,6*L,6)
yz_d = np.arange(3,6*L,6)
yz = np.sort(np.concatenate([yz_u,yz_d]))          
xy_u = np.arange(4,6*L,6)
xy_d = np.arange(5,6*L,6)
xy = np.sort(np.concatenate([xy_u,xy_d]))
#---------------------------------------        

kBT_ext = 3 # How many kB*T intervals to take above the Fermi level
fLim = 1.01 # this factor increases a bit the value of HOL to set the lower bound in eigh method
for hk in Hk: 
    HT = hk + Hv
    # It tries to use scipy routine, which is faster. Otherwise, numpy routine will be used.
    try:
        w, v = LA.eigh(HT, lower = True, overwrite_a = True, subset_by_value = (HOL*fLim, ef + kBT_ext*BPM.eVtoJ*BPM.kB*T))
    except:
        w0, v = np.linalg.eigh(HT, UPLO = 'L')
        w = w0[np.where(w0<(ef + kBT_ext*BPM.eVtoJ*BPM.kB*T))]
    Laux = []
    Oc_xy, Oc_yz, Oc_zx = [], [], []
    for i in range(len(w)):
        Laux.append(BPM.Quasi2DHamiltonian.FermiOcc(w[i])*BPM.Quasi2DHamiltonian.SumInCell(v[:, i],V.L))
        Oc_xy.append(BPM.Quasi2DHamiltonian.FermiOcc(w[i])*BPM.Quasi2DHamiltonian.SumInCell(v[:, i][xy],V.L))
        Oc_yz.append(BPM.Quasi2DHamiltonian.FermiOcc(w[i])*BPM.Quasi2DHamiltonian.SumInCell(v[:, i][yz],V.L))
        Oc_zx.append(BPM.Quasi2DHamiltonian.FermiOcc(w[i])*BPM.Quasi2DHamiltonian.SumInCell(v[:, i][zx],V.L))

    Lrho.append(np.sum(np.array(Laux, dtype = object), axis = 0))
    Lrho_xy.append(np.sum(np.array(Oc_xy, dtype = object), axis = 0))
    Lrho_yz.append(np.sum(np.array(Oc_yz, dtype = object), axis = 0))
    Lrho_zx.append(np.sum(np.array(Oc_zx, dtype = object), axis = 0))
    Oc_xy.clear()
    Oc_yz.clear()
    Oc_zx.clear()
    Laux.clear()

rho = np.sum(np.array(Lrho, dtype = object), axis = 0)/(Nk*Nk) #charge density in e/uc
rho_xy = np.sum(np.array(Lrho_xy, dtype = object), axis = 0)/(Nk*Nk)# xy charge density in e/uc
rho_yz = np.sum(np.array(Lrho_yz, dtype = object), axis = 0)/(Nk*Nk)# yz charge density in e/uc
rho_zx = np.sum(np.array(Lrho_zx, dtype = object), axis = 0)/(Nk*Nk)# zx charge density in e/uc

printlog('Done!')
printlog('Time spent: ' + "{:.2f}".format(t.time()-t2) + ' seg')

printlog('\n')
# plotting the result
PartialChargePlotter(rho, rho_xy, rho_yz, rho_zx, data['DENSITY_ANALYSIS'])
printlog('\n')

end = t.time()
printlog('Total time spent: ' + "{:.2f}".format(end-start) + '  seg')
printlog('\n')
printlog('============================================================================')
now2 = dt.datetime.now()
printlog('Finishing on ' + now2.strftime("%d%b%Y") + ' at ' + now2.strftime("%H:%M:%S"))
printlog('CALCULATION DONE!')
printlog('============================================================================')
printlog('\n') 
plt.show()

