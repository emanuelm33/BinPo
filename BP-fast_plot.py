#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 BinPo, TB-Poisson Solver for 2DES
 Copyright (C) 2021 BinPo Team
 
 BP-fast_plot.py is part of BinPo.
 
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
     FAST PLOTTER TOOL FOR THE SC SOLUTION
     BP-fast_plot.py component allows for quickly plotting the 
     SC potential and the electron density.
"""

__author__ = "Emanuel A. Martinez"
__email__ = "emanuelm@ucm.es"
__copyright__ = "Copyright (C) 2021 BinPo Team"
__version__ = 1.0
__date__ = "January 14, 2022"

import numpy as np
import matplotlib.pyplot as plt
import argparse
plt.style.use('ggplot')

description = 'FAST PLOTTER TOOL FOR THE SC SOLUTION'

parser = argparse.ArgumentParser()
req_name = parser.add_argument_group('required argument')
req_name.add_argument('-id', type = str, help = 'Identifier to plot V and electron density.', required = True)
args = parser.parse_args()

identifier = args.id
#--------------------------------------------------------------------------------------------------------------------------
# MAIN
#####################################################################################
print('\n')
print('=========================================================================')
print('=========================================================================')
print('                                BinPo                                    ')
print('=========================================================================')
print('=========================================================================')
print('\tWELCOME TO THE TIGHT-BINDING POISSON CODE FOR 2DES!')
print('\n')
print('Author: Emanuel A. Martinez, Universidad Complutense de Madrid, Spain.')
print('\n')
print('---------------------------------------------------------------------')
print('\t\tFAST PLOTTER TOOL FOR THE SC SOLUTION')
print('---------------------------------------------------------------------')
print('\n')
print('BinPo, TB-Poisson Solver for 2DES\n'\
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
print('\n')


TBP = np.loadtxt(identifier + '/' + identifier + '_SCP.dat')

fig, (ax1,ax2) = plt.subplots(nrows= 2, ncols = 1, figsize = (6,6.5), sharex = 'all')

ax1.set_ylabel('electrons/plane', size = 15)
ax2.set_ylabel('V [eV]', size = 15)
ax2.set_xlabel('planes', size = 15)

ax1.set_xlim(-2,len(TBP.T[0])+1)

plt.subplots_adjust(left = 0.17, hspace = 0.15, bottom = 0.1, top = 0.95)

d = np.zeros_like(TBP.T[0])
ax1.plot(TBP.T[0],TBP.T[2], lw = 1.6, c = 'teal', marker = 'o', ms = 5, zorder = 10)
ax1.fill_between(TBP.T[0], TBP.T[2], where = TBP.T[2] >= d, interpolate = True, color= 'teal', alpha = 0.3)
ax2.plot(TBP.T[0],TBP.T[1], lw = 1.6, c = 'darkred', marker = 'o', ms = 5, zorder = 10)
print('---------------------------------------------------------------------')
print('\n')
print("-------------------------------------------")
print("Total free charge: " + "{:.4f}".format(np.sum(TBP.T[2])) + '*|e|')
print("-------------------------------------------")

plt.savefig(identifier + '/' + identifier + '_fplot.png', dpi = 300)
print('============================================================================')
print('PLOT DONE!')
print('============================================================================')
print('\n') 
plt.show()

