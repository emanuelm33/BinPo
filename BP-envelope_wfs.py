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
     wavefunctions around the gamma point.
"""

__author__ = "Emanuel A. Martinez"
__email__ = "emanuelm@ucm.es"
__copyright__ = "Copyright (C) 2021 BinPo Team"
__version__ = 1.0
__date__ = "January 14, 2022"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pth
import BPmodule as BPM
import time as t
import datetime as dt
import argparse
import yaml
import logging

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
bc1 = params['BC1_topmost_layer']
bc2 = params['BC2_in_bulk']
ef = params['Fermi_level'] #Fermi level in eV as LUL + dE
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
    # Unfolding dictionary data4
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
    save_plot = data4['SAVING']['save_data']
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
        
    ax.axhline(0,-2,L+1,color = 'r', ls = '-.', zorder = 250)
    ax.plot(z_pl, pot_z, lw = lw_v, marker = 'o', ms = ms, color = c_v, zorder = 230, label = 'SC-potential')
    
    for n in range(N):
        ax.plot(wfs_array[0], wfs_array[n+1], lw = lw_wfs, zorder = 200-2*n, alpha = ALPHA, color = c_wfs)

    if use_legends == True:
        handles = [pth.Circle((0,10), 0.1, color = c_v), pth.Circle((0,10), 0.1, color = c_wfs)]
        ax.legend(handles, ['SC-potential','envelope wfs'], loc = 'best', framealpha = 0.6, facecolor = 'w', fontsize = legend_size)
       
    if save_data == True:
        printlog('Saving data...')
        np.savetxt(identifier + '/' + identifier + '_env-wfs.dat', wfs_array)
    if save_plot == True:
        printlog('Saving plot...')
        plt.savefig(identifier + '/' + identifier + '_env-wfs' + pformat, dpi = resol)
    return None
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
printlog('\t   ENVELOPE WAVEFUNCTIONS AT GAMMA POINT')
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
printlog('\n')
#-------------------------------------------------------------------------------------
# This defines a very little area around the gamma point, with a sampling of 5x5 points
# You can modify these values if necessary
Kmesh = BPM.Kmeshgrid(5, scale = 0.001) 
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

printlog("Transforming <0w|H|Rw'> to k-space...")
t0 = t.time()

T_7 = BPM.Hopping2D(Kmesh,D_7)
T_6 = BPM.Hopping2D(Kmesh,D_6)
T_5 = BPM.Hopping2D(Kmesh,D_5)
T_4 = BPM.Hopping2D(Kmesh,D_4)
T_3 = BPM.Hopping2D(Kmesh,D_3)
T_2 = BPM.Hopping2D(Kmesh,D_2)
T_1 = BPM.Hopping2D(Kmesh,D_1)
T0 = BPM.Hopping2D(Kmesh,D0)

BPM.Quasi2DHamiltonian.set_parameters(T, ef, len(Kmesh))
H = BPM.Quasi2DHamiltonian(T_7,T_6, T_5, T_4, T_3, T_2, T_1, T0, L)
HK = H.HamiltonianTensor()
V = TBP.T[1]


printlog('Initializing potential energy...')
V = BPM.PotentialEnergy(L, bc1, bc2)
V.update_potential(TBP.T[1])
Hv = V.to_tensor()
printlog('Done!')
printlog('\n')

Lenerg, Laux = [], []
env_dict = {}
for n in range(N):
    env_dict[str(n)] = np.zeros(L+1)

for n in range(N):
    for hk in HK: 
        HT = hk + Hv
        w, v = np.linalg.eigh(HT, UPLO = 'L')
        Laux.append(w[n])
        aux = BPM.Quasi2DHamiltonian.SumInCell(v[:, n], L)/len(HK)
        env_dict[str(n)] += np.insert(aux,0,0)
    Lenerg.append(np.sum(np.array(Laux))/len(HK))
    Laux.clear()
    printlog('Wavefunction number ' + str(n + 1) + ' : done!')

Lout = np.zeros((N+1,L+1))
Lout[0] = np.arange(-1,L,1)
for n in range(N):
    Lout[n+1] = env_dict[str(n)]*I_f + (Lenerg[n] - ef)
    
printlog('\n')
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
