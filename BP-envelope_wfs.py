#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 BinPo, TB-Poisson Solver for 2DES
 Copyright (C) 2021 BinPo Team
 
 BP-envelope_wfs.py is part of BinPo.
 
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
     ENVELOPE WAVEFUNCTIONS CALCULATION AT GAMMA POINT
     BP-envelope_wfs.py allows for obtaining and visualizing the envelope
     wavefunctions around the gamma point. At present, this component is 
     only available for cubic systems with t2g manifold, such as STO and KTO.
"""

__author__ = "Emanuel A. Martinez"
__email__ = "emanuelm@ucm.es"
__copyright__ = "Copyright (C) 2021 BinPo Team"
__version__ = 1.1
__date__ = "August 9, 2022"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pth
import BPmodule as BPM
import time as t
import datetime as dt
import argparse
import yaml
import logging
import os 

# Loading the configuration files to set parameters not defined by terminal
#--------------------------------------------------------------------------------------------------------------------------
with open('./config_files/envelope_wfs.yaml', 'r') as f:
    try:
        data = yaml.load(f, Loader=yaml.FullLoader)
    except yaml.YAMLError as exc:
        print ("Error in configuration file: ", exc)
#--------------------------------------------------------------------------------------------------------------------------
# Messages to print by terminal when checking the description with -h or --help command
msj_id = 'Identifier to compute the envelope wavefunctions at gamma point. String.'
msj_sf = 'Scale factor to regulate the intensity of the envelope wavefunctions. Float.'
msj_nw = 'Number of envelope wavefunctions to be computed. Integer.'
msj_xy = 'xy limits for the plot. Float (list). Introduce it as x_min x_max y_min y_max.'
msj_aa = 'Opacity of the envelope wfs. curves. Float. It must be between 0.0 (transparent) and 1.0 (solid).'
description = "ENVELOPE WAVEFUNCTIONS CALCULATION AT GAMMA POINT."
epilog = "By default, arguments not specified are taken from './config_files/envelope_wfs.yaml'."
#--------------------------------------------------------------------------------------------------------------------------
# Passing arguments by terminal and setting their default values
parser = argparse.ArgumentParser(description = description, epilog = epilog)
parser.add_argument('-id', type = str, default = data['ENVELOPE_WAVEFUNCTIONS']['identifier'], help = msj_id)
parser.add_argument('-sf', type = float, default = data['ENVELOPE_WAVEFUNCTIONS']['intensity_factor'], help = msj_sf)
parser.add_argument('-nw', type = int, default = data['ENVELOPE_WAVEFUNCTIONS']['number_of_wavefunctions'], help = msj_nw)
parser.add_argument('-aa', type = float, default = data['ENVELOPE_WAVEFUNCTIONS']['PLOT_ADJUST']['alpha'], help = msj_aa)
parser.add_argument('-xy', nargs = 4, type = float, help = msj_xy)
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
I_f = args.sf
N = args.nw
ALPHA = args.aa
#--------------------------------------------------------------------------------------------------------------------------
# Loading more parameters from the files
L = params['number_of_planes'] # number of planes
material = params['material'] # material
face = params['crystal_face'] # crystal face
a0 = params['lattice_parameter'] # lattice parameter of cubic structure 
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
logger1 = logging.getLogger("BP-env_wfs")
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
def EnvelopePlot(z_pl, pot_z, wfs_array, data4):
    # unfolding data4 dictionary to load several details
    # for the matplotlib plot
    plotstyle = data4['PLOT_ADJUST']['plotstyle']
    fig_size = data4['PLOT_ADJUST']['fig_size']
    c_v = data4['PLOT_ADJUST']['color_V']
    c_wfs = data4['PLOT_ADJUST']['color_wfs']
    ms = data4['PLOT_ADJUST']['markersize']
    l4, b4, r4, t4 = data4['PLOT_ADJUST']['axis_adjust']
    ptitle = data4['PLOT_ADJUST']['title']
    ptitle_size = data4['PLOT_ADJUST']['title_size']
    lw_wfs = data['ENVELOPE_WAVEFUNCTIONS']['PLOT_ADJUST']['linewidth_wfs']
    lw_v = data['ENVELOPE_WAVEFUNCTIONS']['PLOT_ADJUST']['linewidth_V']
    
    xlabel = data4['LABELS']['xlabel']
    xfontsize = data4['LABELS']['xfontsize']
    ylabel = data4['LABELS']['ylabel']
    yfontsize = data4['LABELS']['yfontsize']
    tick_size = data4['LABELS']['ticksize']
    use_legends = data4['LABELS']['legends']
    legend_size = data4['LABELS']['legend_size']
    
    save_data = data4['SAVING']['save_data']
    save_plot = data4['SAVING']['save_plot']
    pformat = data4['SAVING']['format']
    resol = data4['SAVING']['dpi']
    
    plt.style.use(plotstyle)    
    fig, ax = plt.subplots(figsize = fig_size)
    if args.xy:
        ax.set_xlim(args.xy[0], args.xy[1])
        ax.set_ylim(args.xy[2], args.xy[3])
    ax.set_xlabel(xlabel, fontsize = xfontsize)
    ax.set_ylabel(ylabel, fontsize = yfontsize)
    ax.tick_params(labelsize = tick_size)
    if len(ptitle) > 0:
        ax.set_title(ptitle, fontsize = ptitle_size)
    
    plt.subplots_adjust(left = l4, bottom = b4, right = r4, top = t4)    
    ax.axhline(0,-2,L+1,color = 'r', ls = '-.', zorder = -40)
    ax.plot(z_pl, pot_z, lw = lw_v, marker = 'o', ms = ms, color = c_v, zorder = -50, label = 'SC-potential')
    
    for n in range(N):
        ax.plot(wfs_array[0], wfs_array[n+1], lw = lw_wfs, zorder = -5000+2*n, alpha = ALPHA, color = c_wfs)

    if use_legends == True:
        handles = [pth.Circle((0,10), 0.1, color = c_v), pth.Circle((0,10), 0.1, color = c_wfs)]
        ax.legend(handles, ['SC-potential','envelope wfs'], loc = 'best', framealpha = 0.8, facecolor = 'w', fontsize = legend_size)
        
    if save_data == True:
        printlog('Saving data...')
        np.savetxt(identifier + '/' + identifier + '_env-wfs.dat', wfs_array)
    if save_plot == True:
        printlog('Saving plot...')
        plt.savefig(identifier + '/' + identifier + '_env-wfs' + pformat, dpi = resol)
    return None
#--------------------------------------------------------------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------------------------
# NOTE : At moment, this component is only available for cubic ABO3 systems
# with t2g manifold, such as STO and KTO   
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
printlog('\t   ENVELOPE WAVEFUNCTIONS AT GAMMA POINT')
printlog('---------------------------------------------------------------------')
printlog('\n')
printlog('BinPo, TB-Poisson Solver for 2DES\n'\
 'Copyright (C) 2021 BinPo Team\n\n'\
 'BP-envelope_wfs.py is part of BinPo.\n\n'\
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
#-------------------------------------------------------------------------------------
# This defines a very little area around the gamma point, with a sampling of samp_pts x samp_pts 
# points. You can modify these values if necessary
BZ1_frac = 1e-3
samp_pts = 5
Kmesh = BPM.Kmeshgrid(samp_pts, scale = BZ1_frac) 
#------------------------------------------------------------------------------------- 
now = dt.datetime.now()
printlog('DETAILS:')
printlog('\tIdentifier: ' + identifier)
printlog('\tSurface : ' + material + '(' + face + ')')
printlog('\tNumber of planes: ' + str(L))
printlog('\tTemperature: ' + str(T) + ' K')
printlog('\tFermi level: ' + "{:.5f}".format(ef) + ' eV')
printlog('\tNumber of wavefunctions: ' + str(N))
printlog('\n')
printlog('---------------------------------------------------------------------')
printlog('Starting on ' + now.strftime("%d%b%Y") + ' at ' + now.strftime("%H:%M:%S"))
printlog('---------------------------------------------------------------------')
printlog('\n')

if manifold != 't2g':
    printlog('\n')
    printlog('Orbital decomposition of the electron density onto an arbitrary set of MLWFs is not implemented yet.', 'e')
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
BPM.Quasi2DHamiltonian.set_parameters(T, ef, len(Kmesh), HOL)
printlog('Done!')
printlog('Time spent: ' + "{:.2f}".format(t.time()-t0) + ' seg')
printlog('\n')
#############################################################################################################
HK = H.HamiltonianTensor(NWANN) # creating Hamiltonian tensor

printlog('Initializing potential energy...')
V = BPM.PotentialEnergy(L, bc1, bc2) # initializing potential instance
V.update_potential(TBP.T[1]) # updating the potential instance with the SC potential values
Hv = V.to_tensor() # creating potential energy tensor
printlog('Done!')
printlog('\n')

Lenerg, Laux = [], [] # auxiliary lists to save the eigenvalues and the envelope wfcs
env_dict = {} # dictionary to cumulate the envelope wfcs
for n in range(N):
    env_dict[str(n)] = np.zeros(L+1)

for n in range(N): # getting the Schrodinger equation for the N envelope wfcs required
    for hk in HK: 
        HT = hk + Hv
        w, v = np.linalg.eigh(HT, UPLO = 'L')
        Laux.append(w[n])
        aux = BPM.Quasi2DHamiltonian.SumInCell(v[:, n], L)/len(HK)
        env_dict[str(n)] += np.insert(aux,0,0)
    Lenerg.append(np.sum(np.array(Laux))/len(HK))
    Laux.clear()
    printlog('Wavefunction number ' + str(n + 1) + ' : done!')

# preparing the envelope wfcs outputs, by refering them to ef and regulating the size 
# with the intensity factor I_f
Lout = np.zeros((N+1,L+1))
Lout[0] = np.arange(-1,L,1)
for n in range(N):
    Lout[n+1] = env_dict[str(n)]*I_f + (Lenerg[n] - ef)
    
printlog('\n')
# plotting the solution along with the SC-potential to visualize the quantum-well
EnvelopePlot(TBP.T[0], TBP.T[1], Lout, data['ENVELOPE_WAVEFUNCTIONS'])
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
