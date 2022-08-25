#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 BinPo, TB-Poisson Solver for 2DES
 Copyright (C) 2021 BinPo Team
 
 BP-bands.py is part of BinPo.
 
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
     BANDSTRUCTURE CALCULATION
     BP-bands.py component computes the band structure for a 2DES along a
     selected path between high symmetry points in the IBZ1.
"""

__author__ = "Emanuel A. Martinez"
__email__ = "emanuelm@ucm.es"
__copyright__ = "Copyright (C) 2021 BinPo Team"
__version__ = 1.0
__date__ = "January 20, 2022"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pth
import BPmodule as BPM
import time as t
import datetime as dt
import argparse
from matplotlib.colors import Normalize
import yaml
import logging

# Loading the configuration files to set parameters not defined by terminal
#--------------------------------------------------------------------------------------------------------------------------
with open('./config_files/bands.yaml', 'r') as f:
    try:
        data = yaml.load(f, Loader=yaml.FullLoader)
    except yaml.YAMLError as exc:
        print("Error in configuration file: ", exc)
#--------------------------------------------------------------------------------------------------------------------------
# Messages to print by terminal when checking the description with -h or --help command
msj_f = 'Identifier to compute the bandstructure.'
msj_path = 'Path in the IBZ1 to perform the bandstructure calculation. String. Special points for (100) direction are G,X and M, for (110)'\
' direction are G,X,Y,S and for (111) direction are G,K,M.'
msj_nb = "Number of bands to compute. Integer. It must be between 1 and 6*number_of_planes."
msj_kp = 'Number of k-points along the path. Integer.'
msj_tk = "Task to perform. Integer. Options are 0:'total_bandstructure', 1:'orbital_character', 2:'plane_projection' and 3:'all'."
msj_pi = "Initial plane to be considered in planes projection. Integer. Allowed values are between 0 and L-1. Used only if tk = 2 or 3."
msj_pf = "Final plane to be considered in planes projection. Integer. Allowed values are between 'pi' and L. Used only if tk = 2 or 3."
msj_xy = 'Limits for the plot. Float(list). Write it as x_min x_max y_min y_max. To modify separately each type of plot use bands.yaml file.'
description = "BANDSTRUCTURE CALCULATION."
epilog = "By default, arguments not specified are taken from './config_files/bands.yaml'."
#--------------------------------------------------------------------------------------------------------------------------
# Passing arguments by terminal and setting their default values
parser = argparse.ArgumentParser(description = description, epilog = epilog)
parser.add_argument('-id', type = str, default = data['BAND_STRUCTURE']['identifier'], help = msj_f)
parser.add_argument('-ph', type = str, default = data['BAND_STRUCTURE']['path'], help = msj_path)
parser.add_argument('-nb', type = int, default = data['BAND_STRUCTURE']['num_bands'], help = msj_nb)
parser.add_argument('-kp', type = int, default = data['BAND_STRUCTURE']['number_of_kpoints'], help = msj_kp)
parser.add_argument('-tk', type = int, default = data['BAND_STRUCTURE']['bands_task'], help = msj_tk)
parser.add_argument('-pi', type = int, default = data['BAND_STRUCTURE']['initial_plane'], help = msj_pi)
parser.add_argument('-pf', type = int, default = data['BAND_STRUCTURE']['final_plane'], help = msj_pf)
parser.add_argument('-xy', nargs = 4, type = float, default = data['BAND_STRUCTURE']['TOTAL_BANDS']['PLOT_ADJUST']['xy_limits'], help = msj_xy)
args=parser.parse_args()
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
band_task = args.tk
k_points = args.kp
str_path = (args.ph).upper() # uppercase the path
N_BANDS = args.nb
plane1 = args.pi # initial plane considered in projection
plane2 = args.pf # final plane considered in projection

# if xy_lim for total bandstructure is introduced by terminal, automatically the xy_lim for the other two plots are copied 
XY_LIM1 = np.fromiter(args.xy, dtype = np.float64)
if args.xy:
    XY_LIM2 = XY_LIM1
    XY_LIM3 = XY_LIM1
else:
    XY_LIM2 = np.fromiter(data['BAND_STRUCTURE']['ORBITAL_CHARACTER']['xy_limits'], dtype = np.float64) #orbital_character
    XY_LIM3 = np.fromiter(data['BAND_STRUCTURE']['PLANE_PROJECTION']['xy_limits'], dtype = np.float64) #plane_projection
#--------------------------------------------------------------------------------------------------------------------------
# Loading more parameters from the files
L = params['number_of_planes'] # number of planes
material = params['material']
face = params['crystal_face'] # crystal face
a0 = params['lattice_parameter'] # lattice parameter of cubic structure 
T = params['temperature'] # temperature in K
bc1 = params['BC1_topmost_layer']
bc2 = params['BC2_in_bulk']
method = data['BAND_STRUCTURE']['Total_Hk_method'] # method to create the Hamiltonian tensor
F_LEVEL = params['Fermi_level'] #Fermi level in eV as LUL + dE
REF_KPOINT = data['BAND_STRUCTURE']['reference_kpoint']
#--------------------------------------------------------------------------------------------------------------------------
# Creating a logger object for the current program
log_path = identifier + '/' + identifier + '.log'
logger1 = logging.getLogger("BP-bands")
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
# PLOTTING FUNCTIONS DEFINITION
#--------------------------------------------------------------------------------------------------------------------------
# simple band plotter
def Band_Plotter(bandout, pathAxis, data1):
    # Unfolding dictionary data1
    plotstyle = data1['PLOT_ADJUST']['plotstyle']
    lc = data1['PLOT_ADJUST']['linecolor']
    lw = data1['PLOT_ADJUST']['linewidth']
    fig_size = data1['PLOT_ADJUST']['fig_size']
    l1, b1, r1, t1 = data1['PLOT_ADJUST']['axis_adjust']
    ptitle = data1['PLOT_ADJUST']['title']
    ptitle_size = data1['PLOT_ADJUST']['title_size']
    shadow = data1['PLOT_ADJUST']['shadow_above_Ef']
    shadow_color= data1['PLOT_ADJUST']['shadow_color']
        
    xlabel = data1['LABELS']['xlabel']
    xfontsize = data1['LABELS']['xfontsize']
    ylabel = data1['LABELS']['ylabel']
    yfontsize = data1['LABELS']['yfontsize']
    tick_size = data1['LABELS']['ticksize']

    save_data = data1['SAVING']['save_bands']
    save_plot = data1['SAVING']['save_plot']
    pformat = data1['SAVING']['format']
    resol = data1['SAVING']['dpi']

    plt.style.use(plotstyle)
    fig1, ax = plt.subplots(figsize = fig_size)
    index = list(np.where(np.array(pathAxis[2]) == REF_KPOINT)[0])
    xshift = pathAxis[1][index[0]]
    Lout = []
    Lout.append(pathAxis[0]-xshift)
    for i in range(N_BANDS):     
        ax.plot(pathAxis[0]-xshift, bandout.T[i] - F_LEVEL, lw = lw, c = lc, zorder = -30)
        if save_data == True:
            Lout.append(bandout.T[i] - F_LEVEL)
    ax.set_xlabel(xlabel, size = xfontsize)
    ax.set_ylabel(ylabel, size = yfontsize)
    ax.set_xlim(XY_LIM1[0], XY_LIM1[1])
    ax.set_ylim(XY_LIM1[2], XY_LIM1[3])
    ax.tick_params(labelsize = tick_size)
    plt.subplots_adjust(left = l1, bottom = b1, right = r1, top = t1)
    ax.hlines(0,-10,10, color = 'k', lw = 2, ls = '-.', zorder = -10)
    if len(ptitle) > 0:
        ax.set_title(ptitle, fontsize = ptitle_size)
    
    if shadow != 0.0:
        rec = pth.Rectangle((XY_LIM1[0], 0.0), np.abs(XY_LIM1[1]-XY_LIM1[0]), XY_LIM1[3], color = shadow_color, alpha = shadow, zorder = -20)
        ax.add_patch(rec)
        
    if save_data == True:
        printlog('Saving file...')
        np.savetxt(identifier + '/' + identifier + '_bands_' + str_path + '.dat', np.array(Lout).T)
    if save_plot == True:
        printlog('Saving plot...')
        plt.savefig(identifier + '/' + identifier + '_bands_' + str_path + pformat, dpi = resol)
    fig1.show()  

#--------------------------------------------------------------------------------------------------------------------------
# band structure including orbital character
def Band_Plot_TriColor(TBP, TBPc1, TBPc2, TBPc3, pathAxis, data2):
    # Unfolding dictionary data2
    plotstyle = data2['PLOT_ADJUST']['plotstyle']
    c1, c2, c3 = data2['PLOT_ADJUST']['color_seq'].split(',')
    ps = data2['PLOT_ADJUST']['point_size']
    fig_size = data2['PLOT_ADJUST']['fig_size']
    l2, b2, r2, t2 = data2['PLOT_ADJUST']['axis_adjust']
    ptitle = data2['PLOT_ADJUST']['title']
    ptitle_size = data2['PLOT_ADJUST']['title_size']
    shadow = data2['PLOT_ADJUST']['shadow_above_Ef']
    shadow_color = data2['PLOT_ADJUST']['shadow_color']
    
    relsize = data2['COLOR_TRIANGLE']['proportion']
    loc = data2['COLOR_TRIANGLE']['location']
    pad = data2['COLOR_TRIANGLE']['padding']
    font_size = data2['COLOR_TRIANGLE']['fontsize']
    
    xlabel = data2['LABELS']['xlabel']
    xfontsize = data2['LABELS']['xfontsize']
    ylabel = data2['LABELS']['ylabel']
    yfontsize = data2['LABELS']['yfontsize']
    tick_size = data2['LABELS']['ticksize']

    save_data = data2['SAVING']['save_bands']
    save_plot = data2['SAVING']['save_plot']
    pformat = data2['SAVING']['format']
    resol = data2['SAVING']['dpi']

    plt.style.use(plotstyle)
    fig2, ax = plt.subplots(figsize = fig_size)
    index = list(np.where(np.array(pathAxis[2]) == REF_KPOINT)[0])
    xshift = pathAxis[1][index[0]]
    Lout = []
    Lout.append(pathAxis[0]-xshift)
    for j in range(N_BANDS):
        colors = np.array([TBPc1.T[j],TBPc2.T[j],TBPc3.T[j]]).T.reshape(-1,3)
        colors_list = BPM.Gen_Colors(colors, c1, c2, c3, Norm = True)
        ax.scatter(pathAxis[0]-xshift, TBP.T[j]-F_LEVEL, marker = 'o', s = ps, color = colors_list, zorder = -30)
        if save_data == True:
            Lout.append(TBP.T[j]-F_LEVEL)
            Lout.append(TBPc1.T[j])
            Lout.append(TBPc2.T[j])
            Lout.append(TBPc3.T[j])
    ax.set_xlabel(xlabel, size = xfontsize)
    ax.set_ylabel(ylabel, size = yfontsize)
    ax.set_xlim(XY_LIM2[0], XY_LIM2[1])
    ax.set_ylim(XY_LIM2[2], XY_LIM2[3])
    ax.tick_params(labelsize = tick_size)
    plt.subplots_adjust(left = l2, bottom = b2, right = r2, top = t2)
    ax.hlines(0,-10,10, color = 'k', ls = '-.', zorder = -10)
    if len(ptitle) > 0:
        ax.set_title(ptitle, fontsize = ptitle_size)
    #legends
    try:
        BPM.Maxwell_triangle(c1, c2, c3, fig_ax = ax, rel_size = relsize, location = loc, border = pad, fontsize = font_size)
    except:
        handles = [pth.Circle((0,10), 0.1, color = c1), pth.Circle((0,10), 0.1, color = c2), pth.Circle((0,10), 0.1, color = c3)]
        ax.legend(handles, ['xy','yz','zx'], loc = 'best', framealpha = 0.6, facecolor = 'w', fontsize = 12)
        printlog('There was a problem with Maxwell triangle to label the colors. Matplotlib patches were used instead!\n', level = 'w')
    
    if shadow != 0.0:
        rec = pth.Rectangle((XY_LIM1[0], 0.0), np.abs(XY_LIM1[1]-XY_LIM1[0]), XY_LIM1[3], color = shadow_color, alpha = shadow, zorder = -20)
        ax.add_patch(rec)
        
    if save_data == True:
        printlog('Saving file...')
        np.savetxt(identifier + '/' + identifier + '_bands_orbchar_' + str_path + '.dat', np.array(Lout).T)
    if save_plot == True:
        printlog('Saving plot...')
        plt.savefig(identifier + '/' + identifier + '_bands_orbchar_' + str_path + pformat, dpi = resol)
    fig2.show()        
#--------------------------------------------------------------------------------------------------------------------------
# mapcolor for projections onto the B-planes

def PlaneProjector_Plot(TBP, TBP1, pathAxis, data3):
    # Unfolding dictionary data3
    plotstyle = data3['PLOT_ADJUST']['plotstyle']
    ps = data3['PLOT_ADJUST']['point_size']
    fig_size = data3['PLOT_ADJUST']['fig_size']
    c_map = data3['PLOT_ADJUST']['colormap']
    back_color = data3['PLOT_ADJUST']['background_color']
    l3, b3, r3, t3 = data3['PLOT_ADJUST']['axis_adjust']
    ptitle = data3['PLOT_ADJUST']['title']
    ptitle_size = data3['PLOT_ADJUST']['title_size']
    shadow = data3['PLOT_ADJUST']['shadow_above_Ef']
    shadow_color = data3['PLOT_ADJUST']['shadow_color']
    
    locbar = data3['COLORBAR']['location']
    textbar = data3['COLORBAR']['textbar']
    font_size = data3['COLORBAR']['fontsize']
    font_color = data3['COLORBAR']['fontcolor']
        
    xlabel = data3['LABELS']['xlabel']
    xfontsize = data3['LABELS']['xfontsize']
    ylabel = data3['LABELS']['ylabel']
    yfontsize = data3['LABELS']['yfontsize']
    tick_size = data3['LABELS']['ticksize']

    save_data = data3['SAVING']['save_bands']
    save_plot = data3['SAVING']['save_plot']
    pformat = data3['SAVING']['format']
    resol = data3['SAVING']['dpi']
       
    plt.style.use(plotstyle)
    fig3, ax = plt.subplots(figsize = fig_size)
    index = list(np.where(np.array(pathAxis[2]) == REF_KPOINT)[0])
    xshift = pathAxis[1][index[0]]
    Lout = []
    Lout.append(pathAxis[0]-xshift)
    for i in range(N_BANDS):     
        im = ax.scatter(pathAxis[0]-xshift, TBP.T[i]-F_LEVEL, marker = 'o', s = ps, c = TBP1.T[i], cmap = c_map, norm = Normalize(0,1), zorder = -30)
        if save_data == True:
            Lout.append(TBP.T[i]-F_LEVEL)
            # Lout.append(TBP1.T[i]-F_LEVEL) # bug detected!!
            Lout.append(TBP1.T[i])
    cax = fig3.add_axes(locbar) 
    clb = fig3.colorbar(im, ax=ax, cax = cax, ticks = [0.0,1.0])
    clb.ax.set_yticklabels(textbar, size = font_size, color = font_color, fontweight = 'bold')
    clb.outline.set_edgecolor('w')
    cax.yaxis.tick_left()
    #clb.set_ticks([])
    ax.set_xlabel(xlabel, size = xfontsize)
    ax.set_ylabel(ylabel, size = yfontsize)
    ax.set_xlim(XY_LIM3[0], XY_LIM3[1])
    ax.set_ylim(XY_LIM3[2], XY_LIM3[3])
    ax.tick_params(labelsize = tick_size)
    ax.set_facecolor(back_color)
    plt.subplots_adjust(left = l3, bottom = b3, right = r3, top = t3)
    ax.hlines(0,-10,10, color = 'w', ls = '-.', zorder = -10)
    
    if shadow != 0.0:
        rec = pth.Rectangle((XY_LIM1[0], 0.0), np.abs(XY_LIM1[1]-XY_LIM1[0]), XY_LIM1[3], color = shadow_color, alpha = shadow, zorder = -20)
        ax.add_patch(rec)
        
    if len(ptitle) > 0:
        ax.set_title(ptitle, fontsize = ptitle_size)
    if save_data == True:
        printlog('Saving file...')
        np.savetxt(identifier + '/' + identifier + '_bands_planeproj_' + str_path + '.dat', np.array(Lout).T)
    if save_plot == True:
        printlog('Saving plot...')
        plt.savefig(identifier + '/' + identifier + '_bands_planeproj_' + str_path + pformat, dpi = resol)
    fig3.show()        

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
now = dt.datetime.now()
printlog('---------------------------------------------------------------------')
printlog('\t\tBANDSTRUCTURE CALCULATION')
printlog('---------------------------------------------------------------------')
printlog('\n')
printlog('BinPo, TB-Poisson Solver for 2DES\n'\
 'Copyright (C) 2021 BinPo Team\n\n'\
 'BP-bands.py is part of BinPo.\n\n'\
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
# Checking one by one if the values for the input parameters are correct
if plane1 < 0 or plane1 > L-1: # Checking if the values for planes projection are correct
    printlog('initial_plane must be between 0 and L-1', level = 'e')
if plane2 < plane1 or plane2 > L:
    printlog('final_plane cannot be less than initial_plane and cannot be greater than L', level = 'e')

if face == '100': # Checking if the band path is an accepted one
    try:
        path = BPM.BandPath(str_path, BPM.CellGenerator(material, face, a0), k_points)
    except:
        printlog('%s is an invalid path. Special points are G,X and M.' % str_path, level = 'e')
if face == '110':
    try:
        path = BPM.BandPath(str_path, BPM.CellGenerator(material, face, a0), k_points)
    except:
        printlog('%s is an invalid path. Special points are G,X,Y and S.' % str_path, level = 'e')
if face == '111':
    try:
        path = BPM.BandPath(str_path, BPM.CellGenerator(material, face, a0), k_points)
    except:
        printlog('%s is an invalid path. Special points are G,K and M.' % str_path, level = 'e')

if method != 'vectorized' and method != 'iterable':
    printlog("%s is an invalid method. Available methods are 'vectorized' and 'iterable'." % method, level = 'e')

if band_task not in [0,1,2,3]:
    printlog('%s is an invalid option' % band_task, level = 'e')
    
if N_BANDS > 6*L:
    printlog('Number of bands required exceed the total bands in the problem!', level = 'e')
if N_BANDS < 1:
    printlog('The number of bands must be between 1 and 6*number_of_planes!', level = 'e')

#--------------------------------------------------------------------------------------------------------
# Define the k-points using the path
pathAxis = path.get_linear_kpoint_axis()
kpts = path.kpts
dict_task = {0: 'total_bandstructure', 1:'orbital_character', 2:'plane_projection', 3:'total_bands/orb_character/plane_proj'}
#--------------------------------------------------------------------------------------------------------
printlog('DETAILS:')
printlog('\tIdentifier: ' + identifier)
printlog('\tSurface: ' + material + '(' + face + ')')
printlog('\tNumber of planes: ' + str(L))
printlog('\tPath: ' + str_path)
printlog('\tK-points: ' + str(k_points))
printlog('\tNumber of bands: ' + str(N_BANDS))
printlog('\tTemperature: ' + str(T) + ' K')
printlog('\tFermi level: ' + "{:.5f}".format(F_LEVEL) + ' eV')
printlog('\tTask: ' + dict_task[band_task])
if band_task == 2 or band_task == 3:
    printlog('\tInitial plane for projection: ' + str(plane1))
    printlog('\tFinal plane for projection: ' + str(plane2))
printlog('\n')

start = t.time() # general time reference
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

D_7 = BPM.Wann_Sep(Z_7)   # Separation of the <Pw|H|Pw'> for the plane P  
D_6 = BPM.Wann_Sep(Z_6)   # and the Wannier functions w, w'
D_5 = BPM.Wann_Sep(Z_5)
D_4 = BPM.Wann_Sep(Z_4)
D_3 = BPM.Wann_Sep(Z_3)
D_2 = BPM.Wann_Sep(Z_2)
D_1 = BPM.Wann_Sep(Z_1)
D0 = BPM.Wann_Sep(Z0)
#--------------------------------------------------------------------------------------
# 2D Fourier transform of the <Pw|H|Pw'> elements
printlog("Transforming <0w|H|Rw'> to k-space...")
T_7 = BPM.Hopping2D(kpts,D_7)    
T_6 = BPM.Hopping2D(kpts,D_6)
T_5 = BPM.Hopping2D(kpts,D_5)
T_4 = BPM.Hopping2D(kpts,D_4)
T_3 = BPM.Hopping2D(kpts,D_3)
T_2 = BPM.Hopping2D(kpts,D_2)
T_1 = BPM.Hopping2D(kpts,D_1)
T0 = BPM.Hopping2D(kpts,D0)
printlog('Done!')    

BPM.Quasi2DHamiltonian.set_parameters(T,F_LEVEL,len(kpts))
H = BPM.Quasi2DHamiltonian(T_7, T_6, T_5, T_4, T_3, T_2, T_1, T0, L)

V = BPM.PotentialEnergy(L, bc1, bc2)
V.update_potential(TBP.T[1])


if method == 'vectorized':
    HK = H.HamiltonianTensor()
if method == 'iterable':
    HT = []
    for i in range(k_points):
        HT.append(H.HamiltonianMatrix(i))
    HK = np.array(HT)

   
if band_task == 0:
    printlog('\n')
    printlog('Band structure calculation...') 
    result0 = BPM.BandCalculation(HK, V.to_tensor(), N_BANDS)
    Band_Plotter(result0, pathAxis, data['BAND_STRUCTURE']['TOTAL_BANDS'])
    printlog('Done!')
       
if band_task == 1:
    printlog('\n')
    printlog('Band structure calculation with orbital character...') 
    result = BPM.BandCalc_OrbitalChar(HK, V.to_tensor(), L, N_BANDS)
    Band_Plot_TriColor(result[0], result[1], result[2], result[3], pathAxis, data['BAND_STRUCTURE']['ORBITAL_CHARACTER'])
    printlog('Done!')
           
if band_task == 2:
    printlog('\n')
    printlog('Band structure calculation with projections onto planes...')
    result2 = BPM.PlaneProjector(HK, V.to_tensor(), L, plane1, plane2, N_BANDS)
    PlaneProjector_Plot(result2[0], result2[1], pathAxis, data['BAND_STRUCTURE']['PLANE_PROJECTION'])
    printlog('Done!')

if band_task == 3:
    printlog('\n')
    printlog('Band structure calculation...') 
    result0 = BPM.BandCalculation(HK, V.to_tensor(), N_BANDS)
    Band_Plotter(result0, pathAxis, data['BAND_STRUCTURE']['TOTAL_BANDS'])
    printlog('Done!')
    printlog('\n')
    printlog('Band structure calculation with orbital character...') 
    result = BPM.BandCalc_OrbitalChar(HK, V.to_tensor(), L, N_BANDS)
    Band_Plot_TriColor(result[0], result[1], result[2], result[3], pathAxis, data['BAND_STRUCTURE']['ORBITAL_CHARACTER'])
    printlog('Done!')
    printlog('\n')
    printlog('Band structure calculation with projections onto planes...')
    result2 = BPM.PlaneProjector(HK, V.to_tensor(), L, plane1, plane2, N_BANDS)
    PlaneProjector_Plot(result2[0], result2[1], pathAxis, data['BAND_STRUCTURE']['PLANE_PROJECTION'])
    printlog('Done!')
  
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

