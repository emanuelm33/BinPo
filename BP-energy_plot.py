#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 BinPo, TB-Poisson Solver for 2DES
 Copyright (C) 2021 BinPo Team
 
 BP-energy_plot.py is part of BinPo.
 
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
     ENERGY SLICES PLOTTER
     BP-energy_plot.py plots the output of BP-energy_slices.py component.
"""

import numpy as np
import matplotlib.pyplot as plt
import BPmodule as BPM
import matplotlib.patches as pth
import argparse
import datetime as dt
import time as t
import yaml
import logging

# Loading the configuration files to set parameters not defined by terminal
#--------------------------------------------------------------------------------------------------------------------------
with open('./config_files/energy_plot.yaml', 'r') as f:
    try:
        data = yaml.load(f, Loader=yaml.FullLoader)
    except yaml.YAMLError as exc:
        print ("Error in configuration file: ", exc)
#--------------------------------------------------------------------------------------------------------------------------
# Messages to print by terminal when checking the description with -h or --help command
msj_id = 'Identifier to plot an energy slice. String.'
msj_in = 'Name of the input file for plotting the energy slice. String.'
msj_ew = 'Energy window. Float. The points shown will range from (energy_cut-energy_window) to energy_cut.'
msj_ps = 'Pointsize for scatter plot. Float.'
msj_xy = 'Limits for the plot. Float (list). Write it as x_min x_max y_min y_max.'
description = "ENERGY SLICES PLOTTER"
epilog = "By default, arguments not specified are taken from './config_files/energy_plot.yaml'."
#--------------------------------------------------------------------------------------------------------------------------
# Passing arguments by terminal and setting their default values
parser = argparse.ArgumentParser(description = description, epilog = epilog)
parser.add_argument('-id', type = str, default = data['ENERGY_PLOTTER']['identifier'], help = msj_id)
parser.add_argument('-i', type = str, default = data['ENERGY_PLOTTER']['input_file'], help = msj_in)
parser.add_argument('-ew', type = float, default = data['ENERGY_PLOTTER']['energy_window'], help = msj_ew)
parser.add_argument('-ps', type = float, default = data['ENERGY_PLOTTER']['PLOT_ADJUST']['point_size'], help = msj_ps)
parser.add_argument('-xy', nargs = 4, type = float, default = data['ENERGY_PLOTTER']['PLOT_ADJUST']['xy_limits'], help = msj_xy)
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
#--------------------------------------------------------------------------------------------------------------------------
# Copying the rest of the parameters updatable by terminal
# The uppercase is because of their treatment as global variables here!
XY_LIM = args.xy # xy limits of the plot
WNW = args.ew # energy window
PTS = args.ps # point size
IN_FILE = args.i # name of input file
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
logger1 = logging.getLogger("BP-eplot")
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

# PLOTTING FUNCTION DEFINITION
#--------------------------------------------------------------------------------------------------------------------------

def ESlicesPlot_OrbChar(x, y, z, Cxy, Cyz, Czx, face, data0):
    # Unfolding data0 dictionary
    orbchar = data0['orbital_char']
    c1, c2, c3 = data0['color_seq'].split(',')
    
    plotstyle = data0['PLOT_ADJUST']['plotstyle']
    fig_size = data0['PLOT_ADJUST']['fig_size']
    l0, b0, r0, t0 = data0['PLOT_ADJUST']['axis_adjust']
    ptitle = data0['PLOT_ADJUST']['title']
    ptitle_size = data0['PLOT_ADJUST']['title_size']
    
    xlabel = data0['LABELS']['xlabel']
    xfontsize = data0['LABELS']['xfontsize']
    ylabel = data0['LABELS']['ylabel']
    yfontsize = data0['LABELS']['yfontsize']
    tick_size = data0['LABELS']['ticksize']
    
    relsize = data0['COLOR_TRIANGLE']['proportion']
    loc = data0['COLOR_TRIANGLE']['location']
    pad = data0['COLOR_TRIANGLE']['padding']
    font_size = data0['COLOR_TRIANGLE']['fontsize']
    
    save_plot = data0['SAVING']['save_plot']
    pformat = data0['SAVING']['format']
    resol = data0['SAVING']['dpi']
    
    
    plt.style.use(plotstyle)
    fig2, ax = plt.subplots(figsize = fig_size) 
    ax.set_aspect(1.0)
    plt.subplots_adjust(left = 0.2)
    z_filt = np.where(np.logical_and(z!=-1, z > (np.max(z) - WNW)))
    X = x[z_filt]
    Y = y[z_filt]
    Cxy_f = Cxy[z_filt]
    Cyz_f = Cyz[z_filt]
    Czx_f = Czx[z_filt]
    colors = np.array([Cxy_f/np.max(Cxy_f), Cyz_f/np.max(Cyz_f), Czx_f/np.max(Czx_f)]).T.reshape(-1,3)
    if orbchar == True:
        colors_list = BPM.Gen_Colors(colors, c1, c2, c3)
        #legends
        try:
            BPM.Maxwell_triangle(c1, c2, c3, fig_ax = ax, rel_size = relsize, location = loc, border = pad, fontsize = font_size)    
        except:
            handles = [pth.Circle((0,10), 0.1, color = c1), pth.Circle((0,10), 0.1, color = c2), pth.Circle((0,10), 0.1, color = c3)]
            ax.legend(handles, ['xy','yz','zx'], loc = 'best', framealpha = 0.6, facecolor = 'w', fontsize = 12)
            printlog('There was a problem with Maxwell triangle to label the colors. Matplotlib patches were used instead!\n', level = 'w')
    else:
        colors_list = c1
    
    if face in ['110','101','011']:
        r2 = np.sqrt(2)
        Rt = np.array([[1/r2,-1/r2,0.0],[0.0,0.0,-1.0],[1/r2,1/r2,0.0]])
    
        b1 = np.array([-0.5,0.5,0.0])
        b2 = np.array([0.0,0.0,1.0])

        T_kmesh_aux = []
        for j in range(len(X)):
            T_kmesh_aux.append(np.dot(Rt,X[j]*b1+Y[j]*b2))
        T_kmesh = np.array(T_kmesh_aux)
        
        X = T_kmesh.T[0]
        Y = T_kmesh.T[1]
    
    if face == '111':
        Rt = np.array([[0.707,-0.707,0],[0.408,0.408,-0.817],[0.578,0.578,0.577]])

        b1 = np.array([-1/3,-1/3,2/3])
        b2 = np.array([2/3,-1/3,-1/3])
    
        T_kmesh_aux = []
        for j in range(len(X)):
            T_kmesh_aux.append(np.dot(Rt,X[j]*b1+Y[j]*b2))
        T_kmesh = np.array(T_kmesh_aux)
        
        X = T_kmesh.T[0]
        Y = T_kmesh.T[1]
              
    ax.set_xlim(XY_LIM[0], XY_LIM[1])
    ax.set_ylim(XY_LIM[2], XY_LIM[3])
    ax.set_xlabel(xlabel, size = xfontsize)
    ax.set_ylabel(ylabel, size = yfontsize)
    ax.tick_params(labelsize = tick_size)
    plt.subplots_adjust(left = l0, bottom = b0, right = r0, top = t0)
    if len(ptitle) > 0:
        ax.set_title(ptitle, fontsize = ptitle_size)
    ax.scatter(X, Y, marker = 'o', s = PTS, color = colors_list, zorder = 50)
    if save_plot == True:
        printlog('Saving plot...')
        if IN_FILE == 'default':
            plt.savefig(identifier + '/' + identifier + '_ES_plot' + pformat, dpi = resol)
        else:
            plt.savefig(identifier + '/' + IN_FILE + '_ES_plot' + pformat, dpi = resol)
    return None
#--------------------------------------------------------------------------------------------------------------------------
# MAIN
######################################################################################
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
printlog('\t\tENERGY SLICES PLOTTER')
printlog('---------------------------------------------------------------------')
printlog('\n')
printlog('BinPo, TB-Poisson Solver for 2DES\n'\
 'Copyright (C) 2021 BinPo Team\n\n'\
 'BP-energy_plot.py is part of BinPo.\n\n'\
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
# Loading and unfolding the output of BP-energy_slices.py
if IN_FILE == 'default':
    try:
        F = np.loadtxt(identifier + '/' + identifier + '_ES.dat')
    except:
        printlog('The input file could not be loaded', level = 'e')
else:
    try:
        F = np.loadtxt(identifier + '/' + IN_FILE + '.dat')
    except:
        printlog('The input file could not be loaded', level = 'e')
x, y, z, Cxy, Cyz, Czx = F.T[0]*2*np.pi/a0, F.T[1]*2*np.pi/a0, F.T[2], F.T[3], F.T[4], F.T[5] 
#--------------------------------------------------------------------------------------------------------------------------
now = dt.datetime.now()
printlog('\n')
printlog('---------------------------------------------------------------------')
printlog('Starting on ' + now.strftime("%d%b%Y") + ' at ' + now.strftime("%H:%M:%S"))
printlog('---------------------------------------------------------------------')
printlog('\n')
printlog('DETAILS:')
printlog('\tIdentifier : ' + identifier)
if IN_FILE == 'default':
    printlog('\tInput file : ' + identifier + '_ES.dat')
else:
    printlog('\tInput file : ' + IN_FILE + '.dat')
printlog('\tSurface : ' + material + '(' + face + ')')
printlog('\tNumber of planes : ' + str(L))
printlog('\tTotal k-points : ' + str(len(x)))
printlog('\tTemperature : ' + str(T) + ' K')
printlog('\tFermi level : ' + "{:.5f}".format(ef) + ' eV')
printlog('\tEnergy window : ' + str(WNW) + ' eV')
printlog('\tOrbital character : ' + str(data['ENERGY_PLOTTER']['orbital_char']))
printlog('\n')

start = t.time() # general time reference

printlog('Plotting the energy slice')
printlog('It could take a few seconds...')
printlog('\n')

ESlicesPlot_OrbChar(x, y, z, Cxy, Cyz, Czx, face, data['ENERGY_PLOTTER'])

end = t.time()
printlog('Total time spent: ' + "{:.2f}".format(end-start) + '  seg') 
printlog('\n')
printlog('============================================================================')
now2 = dt.datetime.now()
printlog('Finishing on ' + now2.strftime("%d%b%Y") + ' at ' + now2.strftime("%H:%M:%S"))
printlog('PLOT DONE!')
printlog('============================================================================')
printlog('\n')
plt.show()

