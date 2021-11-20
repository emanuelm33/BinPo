#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 BinPo, TB-Poisson Solver for 2DES
 Copyright (C) 2021 BinPo Team
 
 BP-preproc.py is part of BinPo.
 
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
     File for pre-processing step.
     The Wannier file is taken from 'WannFolder' and a new folder containing
     the elements of r-space Hamiltonian separated by planes is generated.
     In principle, this component should be executed once per crystal face. 
     After running it the user can perform any SC-potential calculation associated
     to the specific Wannier file.
     This component does not log the information and need two mandatory input
     through the terminal: the material and the crystal face.
"""
import os
import numpy as np
import time as t
import datetime as dt
import argparse
import BPdatabase as BPD

#--------------------------------------------------------------------------------------------------------------------------
# Messages to print by terminal when checking the description with -h or --help command
description = "PREPROCESSING STEP. The Wannier file is taken from 'WannFolder' and a new folder containing"\
    " the elements of r-space Hamiltonian separated by planes is generated."
epilog = "In principle, this component should be executed once per crystal face. After running it"\
    " the user can perform any SC-potential calculation associated to the specific Wannier file."
msj_fc = "Crystal face. This argument is mandatory. String. It could be either '100' ('010', '001') or '110' ('011', '101') or '111'."
msj_mat = "Name of the material. This argument is mandatory. String. It could be 'STO' or 'KTO'."
#--------------------------------------------------------------------------------------------------------------------------
# Passing mandatory arguments by terminal
parser = argparse.ArgumentParser(description = description, epilog = epilog)
req_name = parser.add_argument_group('required argument')
req_name.add_argument('-m', type = str, help = msj_mat, required = True)
req_name.add_argument('-fc', type = str, help = msj_fc, required = True)
args = parser.parse_args()
#--------------------------------------------------------------------------------------------------------------------------
# Loading the values
material = args.m
if material != 'KTO' and material != 'STO':
    raise ValueError('%s is an invalid material!')
face = args.fc
#--------------------------------------------------------------------------------------------------------------------------
filename = BPD.materials[material]['Wannier_file'] # open the Wannier file
F001 = np.loadtxt('./WFolder/'+ filename, skiprows = BPD.materials[material]['skip_rows'])
F001c = np.copy(F001) # security copy of the loaded Wannier file
#--------------------------------------------------------------------------------------------------------------------------
def planes_sep(file, value):
    """ Separate the bulk Hamiltonian into different plane contributions 
    according to <0w|H|Rw'>, where R is the perpendicular component to the 
    selected crystal face and w,w' are the Wannier functions indices."""
    OR = np.logical_or(file.T[5] != 0,file.T[6] != 0)
    indices = np.unique(np.array([np.argwhere(np.logical_and(file.T[2] == value, OR)) for i in file.T[2]]))
    return np.array(file[indices])
#--------------------------------------------------------------------------------------------------------------------------
print('\n')
print('=========================================================================')
print('=========================================================================')
print('                                BinPo                                    ')
print('=========================================================================')
print('=========================================================================')
now = dt.datetime.now()
print('\tWELCOME TO THE TIGHT-BINDING POISSON CODE FOR 2DES!')
print('\n')
print('Author: Emanuel A. Martinez, Universidad Complutense de Madrid, Spain.')
print('\n')
print('-------------------------------------------------------------------------')
print('\t\t\tPRE-PROCESSING STEP')
print('-------------------------------------------------------------------------')
print('\n')
print('BinPo, TB-Poisson Solver for 2DES\n'\
 'Copyright (C) 2021 BinPo Team\n\n'\
 'BP-preproc.py is part of BinPo.\n\n'\
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
print('-------------------------------------------------------------------------')
print('Starting on ' + now.strftime("%d%b%Y") + ' at ' + now.strftime("%H:%M:%S"))
print('-------------------------------------------------------------------------')
print('\n')
print('Surface : ' + material + '(' + face + ')')
print('Performing a separation of planes along normal direction.')
print('It could take a few minutes...')

t0 = t.time()

if face in ['100', '010', '001']:
    
    face = '100'
    Z_7 = planes_sep(F001,-7)
    print('Plane -7 done!')
    Z_6 = planes_sep(F001,-6)
    print('Plane -6 done!')
    Z_5 = planes_sep(F001,-5)
    print('Plane -5 done!')
    Z_4 = planes_sep(F001,-4)
    print('Plane -4 done!')
    Z_3 = planes_sep(F001,-3)
    print('Plane -3 done!')
    Z_2 = planes_sep(F001,-2)
    print('Plane -2 done!')
    Z_1 = planes_sep(F001,-1)
    print('Plane -1 done!')
    Z0 = planes_sep(F001,0)
    print('Plane 0 done!')


elif face in ['110', '101', '011']:
    face = '110'
    #basis change
    C_NB = np.array([[-1,1,0],[0,0,1],[1,1,0]]).T#C[N->B] 
    C_BN = np.linalg.inv(C_NB)#C[B->N]

    f110 = []
    for line in F001c:
        line[:3] = np.around(np.dot(C_BN,line[:3]),2)
        f110.append(line)

    F110 = np.array(f110)
    
    Z_7 = planes_sep(F110,-3.50)
    print('Plane -7 done!')
    Z_6 = planes_sep(F110,-3.00)
    print('Plane -6 done!')
    Z_5 = planes_sep(F110,-2.50)
    print('Plane -5 done!')
    Z_4 = planes_sep(F110,-2.00)
    print('Plane -4 done!')
    Z_3 = planes_sep(F110,-1.50)
    print('Plane -3 done!')
    Z_2 = planes_sep(F110,-1.00)
    print('Plane -2 done!')
    Z_1 = planes_sep(F110,-0.50)
    print('Plane -1 done!')
    Z0 = planes_sep(F110,0.00)
    print('Plane 0 done!')

    
elif face == '111':

    #basis change
    C_NB = np.array([[0,-1,1],[1,-1,0],[1,1,1]]).T#C[N->B] 
    C_BN = np.linalg.inv(C_NB)#C[B->N]

    f111 = []
    for line in F001c:
        line[:3] = np.around(np.dot(C_BN,line[:3]),2)
        f111.append(line)

    F111 = np.array(f111)
    
    Z_7 = planes_sep(F111, -2.33)
    print('Plane -7 done!')
    Z_6 = planes_sep(F111, -2.00)
    print('Plane -6 done!')
    Z_5 = planes_sep(F111, -1.67)
    print('Plane -5 done!')
    Z_4 = planes_sep(F111, -1.33)
    print('Plane -4 done!')
    Z_3 = planes_sep(F111, -1.00)
    print('Plane -3 done!')
    Z_2 = planes_sep(F111, -0.67)
    print('Plane -2 done!')
    Z_1 = planes_sep(F111, -0.33)
    print('Plane -1 done!')
    Z0 = planes_sep(F111, 0.00)
    print('Plane 0 done!')

else:
    raise ValueError('%s is an invalid face. Allowed values are 001 (010,100), 110 (101,011) and 111.' % face)
    
os.system('mkdir Hr' + material + face)
print('Saving files to the Hr' + material + face + ' folder...' )
np.savetxt('Hr' + material + face + '/Z_7.dat', Z_7, fmt = ['%.2f','%.2f','%.2f','%.1f','%.1f','%.18e','%.18e'])
np.savetxt('Hr' + material + face + '/Z_6.dat', Z_6, fmt = ['%.2f','%.2f','%.2f','%.1f','%.1f','%.18e','%.18e'])
np.savetxt('Hr' + material + face + '/Z_5.dat', Z_5, fmt = ['%.2f','%.2f','%.2f','%.1f','%.1f','%.18e','%.18e'])
np.savetxt('Hr' + material + face + '/Z_4.dat', Z_4, fmt = ['%.2f','%.2f','%.2f','%.1f','%.1f','%.18e','%.18e'])
np.savetxt('Hr' + material + face + '/Z_3.dat', Z_3, fmt = ['%.2f','%.2f','%.2f','%.1f','%.1f','%.18e','%.18e'])
np.savetxt('Hr' + material + face + '/Z_2.dat', Z_2, fmt = ['%.2f','%.2f','%.2f','%.1f','%.1f','%.18e','%.18e'])
np.savetxt('Hr' + material + face + '/Z_1.dat', Z_1, fmt = ['%.2f','%.2f','%.2f','%.1f','%.1f','%.18e','%.18e'])
np.savetxt('Hr' + material + face + '/Z0.dat', Z0, fmt = ['%.2f','%.2f','%.2f','%.1f','%.1f','%.18e','%.18e'])
print('Done!')
print('Time spent: ' + "{:.2f}".format(t.time()-t0) + '  seg')
print('\n')
print('=========================================================================')
now2 = dt.datetime.now()
print('Finishing on ' + now2.strftime("%d%b%Y") + ' at ' + now2.strftime("%H:%M:%S"))
print('CALCULATION DONE!')
print('=========================================================================')
print('\n')  
    
