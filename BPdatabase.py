#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 BinPo, TB-Poisson Solver for 2DES
 Copyright (C) 2021 BinPo Team
 
 BPdatabase.py is part of BinPo.
 
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
     This module contains a dictionary with the available materials, for the time being.
     You could add more items with different materials or different DFT treatments. 
     The composition of each dictionary as item of materials main dictionary is:
         
         key = 'name of the perovskite complex oxide ABO3'.
         values = { description : the formula and another info.
                           DATA USED IN DFT CALCULATION
                    XC-functional.
                    pseudopotential_type.
                           DATA OBTAINED FROM DFT CALCULATION
                    lattice_parameter (Angs).
                    lowest_unoccupied_level (eV).
                           DATA FOR LOADING THE WANNIER FILE
                    Wannier_file : name of the file.
                    skip_rows : number of lines to skip when reading the 
                    Wannier file.
                 } 
    NOTE : If you want to add a new item, make sure of using differents keys for the same material. For example,
    to add a new item for SrTiO3 with different DFT XC-functional and/or pseudopotential,
    name the key as 'STO1' or something like this. Otherwise you will obtain an error by 
    repeating keys in the main dictionary 'materials'.
"""
__version__ = "1.0"
__author__ = "Emanuel A. Martinez"

materials = { 'STO' : {'description' : 'cubic SrTiO3',
                       'XC-functional' : 'PBE',
                       'pseudopotential_type' : 'PAW',
                       'lattice_parameter' : 3.9425,
                       'lowest_unoccupied_level' : 11.4685,
                       'Wannier_file' : 'STO_hr.dat',
                       'skip_rows' : 0
                       },


              'KTO' : { 'description' : 'cubic KTaO3',
                        'XC-functional' : 'PBE',
                        'pseudopotential_type' : 'PAW',
                        'lattice_parameter' : 4.0184,
                        'lowest_unoccupied_level' : 10.5168,
                        'Wannier_file' : 'KTO_hr.dat',
                        'skip_rows' : 0
                       },

              # Add here a new key : values item
              # 'name' : { 'description' : ,
              #           'XC-functional' : ,
              #           'pseudopotential_type' : ,
              #           'lattice_parameter' : ,
              #           'lowest_unoccupied_level' : ,
              #           'Wannier_file' : ,
              #           'skip_rows' : 
              #          },
             }

###############################################################################################################################

if __name__ == "__main__":
    print('\n')
    print('=========================================================================')
    print('=========================================================================')
    print('                                BinPo                                    ')
    print('=========================================================================')
    print('=========================================================================')
    print('  WELCOME TO THE TIGHT-BINDING POISSON CODE FOR 2DES!')
    print('\n')
    print('Author: Emanuel A. Martinez, Universidad Complutense de Madrid, Spain.')
    print('\n')
    print('BinPo database module. This is not a runnable file.')
    print('\n')
    print('=========================================================================')
    print('\t DOCUMENTATION')
    print('=========================================================================')
    print(help(__name__))    
    print('=========================================================================')
    print('=========================================================================')

