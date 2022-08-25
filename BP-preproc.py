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

__author__ = "Emanuel A. Martinez"
__email__ = "emanuelm@ucm.es"
__copyright__ = "Copyright (C) 2021 BinPo Team"
__version__ = 1.1
__date__ = "August 9, 2022"

import os
import numpy as np
import time as t
import datetime as dt
import argparse
import BPdatabase as BPD

#--------------------------------------------------------------------------------------------------------------------------
# Messages to print by terminal when checking the description with -h or --help command
description = "PREPROCESSING STEP. The Wannier90 file (W90 file) is taken from 'WFolder' and a new folder containing"\
    " the elements of r-space Hamiltonian separated by planes is generated."
epilog = "In principle, this component should be executed once per material + confinement direction. After running it"\
    " the user can perform any SC-potential calculation associated to the W90 file/material/direction combination."
msj_fc = "Confinement direction 'hkl'. String. This argument is mandatory. For cubic systems, it could be either '100' ('010', '001') or"\
    " '110' ('011', '101') or '111'. For hexagonal systems, it must be '100' ('010', '001')."
msj_mat = "Name of the material. String. This argument is mandatory. It must match the corresponding key in materials dictionary "\
    "of BPdatabase.py module."
#--------------------------------------------------------------------------------------------------------------------------
# Passing mandatory arguments by terminal
parser = argparse.ArgumentParser(description = description, epilog = epilog)
req_name = parser.add_argument_group('required argument')
req_name.add_argument('-mt', type = str, help = msj_mat, required = True)
req_name.add_argument('-cfd', type = str, help = msj_fc, required = True)
args = parser.parse_args()
#--------------------------------------------------------------------------------------------------------------------------
# Loading the values
material = (args.mt).upper() # uppercase the material name
face = args.cfd
#--------------------------------------------------------------------------------------------------------------------------
def planes_sep(file, value):
    """ Separate the bulk Hamiltonian into different plane contributions according
    to <0 w|H|Rxy + Rz w'>, where Rz is parallel to the confinement direction and perpendicular
    to the in-plane elements Rxy, while the w,w' are the Wannier functions indices."""
    sep = []
    for line in file:
    	if line[2] == value:
    		sep.append(line)    	
    return np.array(sep)
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
 'preproc.py is part of BinPo.\n\n'\
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
print('It could take a few seconds...')
#--------------------------------------------------------------------------------------------------------------------------
filename = BPD.materials[material]['Wannier_file'] # opening the W90 file
F001 = np.loadtxt('./WFolder/'+ filename, skiprows = BPD.materials[material]['skip_rows'])
#--------------------------------------------------------------------------------------------------------------------------
sgeom = BPD.materials[material]['unit-cell'] # getting the system geometry
#--------------------------------------------------------------------------------------------------------------------------
# Filtering of elements equal to zero in the W90 file
print('Filtering zero elements...')
Fcopy = []
for i in F001:
    if np.logical_or(i[5] != 0.0, i[6] != 0.0) == True:
    	Fcopy.append(i)

F001f = np.array(Fcopy) # Filtered new file
print('Done!')	
F001c = np.copy(F001f) # security copy (useful for the change of basis)
#--------------------------------------------------------------------------------------------------------------------------
t0 = t.time() # General time reference

# First, the system is classified by geometry, then by confinement direction (face)
planes = [] # auxiliary list to save the <0 w|H|Rxy + Rz w'> elements separated by planes

if sgeom == 'cubic':
    
    if face in ['100', '010', '001']:
        face = '100'
        # the dictionary values indicate the Rz position (in lattice fraction) where the separation
        # by planes of the <0 w|H|Rxy + Rz w'> elements must be performed.
        dict100 = {0:-7, 1:-6, 2:-5, 3:-4, 4:-3, 5:-2, 6:-1, 7:0}
        # performing the separation and saving it to the planes list
        for i in dict100.values():
             planes.append(planes_sep(F001f,i))
             print('Plane ' + str(i) + ' done!')
           
    elif face in ['110', '101', '011']:
        face = '110'
        # basis change from (001) to (110) direction
        C_NB = np.array([[-1,1,0],[0,0,1],[1,1,0]]).T#C[N->B] 
        C_BN = np.linalg.inv(C_NB)#C[B->N]
    
        f110 = []
        for line in F001c:
            line[:3] = np.around(np.dot(C_BN,line[:3]),2)
            f110.append(line)
    
        F110 = np.array(f110)
        # the dictionary values indicating the Rz for plane separation
        dict110 = {0:-3.50, 1:-3.00, 2:-2.50, 3:-2.00, 4:-1.50, 5:-1.00, 6:-0.50, 7:0.00}
        
        for key,val in dict110.items():
            planes.append(planes_sep(F110,val))
            print('Plane ' + str(key) + ' done!')
            
    elif face == '111':
        # basis change from (001) to (111) direction
        C_NB = np.array([[0,-1,1],[1,-1,0],[1,1,1]]).T#C[N->B] 
        C_BN = np.linalg.inv(C_NB)#C[B->N]
    
        f111 = []
        for line in F001c:
            line[:3] = np.around(np.dot(C_BN,line[:3]),2)
            f111.append(line)
    
        F111 = np.array(f111)
        # the dictionary values indicating the Rz for plane separation
        dict111 = {0:-2.33, 1:-2.00, 2:-1.67, 3:-1.33, 4:-1.00, 5:-0.67, 6:-0.33, 7:0.00}
        
        for key,val in dict111.items():
            planes.append(planes_sep(F111,val))
            print('Plane ' + str(key) + ' done!')
    
    else:
        raise ValueError('%s is an invalid face. Allowed values are 001 (010,100), 110 (101,011), 111.' % face)


elif sgeom == 'hexagonal':

    if face in ['100', '010', '001']:
        face = '001h' # rename the face for the hexagonal case
        # the dictionary values indicating the Rz for plane separation
        dict001h = {0:-7.00, 1:-6.00, 2:-5.00, 3:-4.00, 4:-3.00, 5:-2.00, 6:-1.00, 7:0.00}
    
        for key,val in dict001h.items():
            planes.append(planes_sep(F001c,val))
            print('Plane ' + str(key) + ' done!')
    else:
        raise ValueError('%s is an invalid face. Allowed values are 001 (010,100).' % face)
    
else:
    raise ValueError("%s is an invalid system geometry. Allowed values are 'cubic' or 'hexagonal'." % sgeom)

# creating a folder to save one-by-one the <0 w|H|Rxy + Rz w'> elements, separated by Rz value from planes list 
os.system('mkdir Hr' + material + face)
print('Saving files to the Hr' + material + face + ' folder...' )

N_planes = 8
for i in range(0, N_planes):
    # if for some Rz there are no elements, an empty file is saved instead
    try:
    	np.savetxt('Hr' + material + face + '/Z' + str(i-(N_planes-1)) + '.dat', planes[i], fmt = ['%.2f','%.2f','%.2f','%d','%d','%.18e','%.18e'])
    except:
    	print('Z' + str(i) + '.dat' + ' file is empty!')

planes.clear() # cleaning the auxiliary list

print('Done!')
print('Time spent: ' + "{:.2f}".format(t.time()-t0) + '  seg')
print('\n')
print('=========================================================================')
now2 = dt.datetime.now()
print('Finishing on ' + now2.strftime("%d%b%Y") + ' at ' + now2.strftime("%H:%M:%S"))
print('CALCULATION DONE!')
print('=========================================================================')
print('\n')  

