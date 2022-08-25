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
     You could add more items for different systems (see ~/BPexamples/Adding_W90files.pdf).
     
"""

__author__ = "Emanuel A. Martinez"
__email__ = "emanuelm@ucm.es"
__copyright__ = "Copyright (C) 2021 BinPo Team"
__version__ = 1.1
__date__ = "August 9, 2022"


materials = { 'STO' : {'description' : 'cubic SrTiO3, Ti t2g manifold',
                       'pseudopotential' : 'PAW PBE full-relativistic from PSlibrary',
                       'unit-cell' : 'cubic',
                       'Wannier_functions_number' : 6,
                       'lattice_parameter' : 3.9425,
                       'highest_occupied_level' : 9.6404,
                       'lowest_unoccupied_level' : 11.4685,
                       'Wannier_file' : 'STO_hr.dat',
                       'skip_rows' : 228,
                       'manifold' : 't2g',
                       },
             
             
             'STOB' : {'description' : 'cubic SrTiO3, Ti t2g and O 2p manifolds',
                       'pseudopotential' : 'PAW PBE full-relativistic from PSlibrary',
                       'unit-cell' : 'cubic',
                       'Wannier_functions_number' : 24,
                       'lattice_parameter' : 3.9425,
                       'highest_occupied_level' : 9.6404,
                       'lowest_unoccupied_level' : 11.4685,
                       'Wannier_file' : 'STOB_hr.dat',
                       'skip_rows' : 228,
                       'manifold' : 'other',
                       },
             
             'STOC' : {'description' : 'cubic SrTiO3, Ti t2g and Ti eg manifolds',
                       'pseudopotential' : 'PAW PBE full-relativistic from PSlibrary',
                       'unit-cell' : 'cubic',
                       'Wannier_functions_number' : 10,
                       'lattice_parameter' : 3.9425,
                       'highest_occupied_level' : 9.6404,
                       'lowest_unoccupied_level' : 11.4685,
                       'Wannier_file' : 'STOC_hr.dat',
                       'skip_rows' : 228,
                       'manifold' : 'other',
                       },
             
              'KTO' : {'description' : 'cubic KTaO3, Ta t2g manifold',
                       'pseudopotential' : 'NC PBE full-relativistic from Pseudo-Dojo',
                       'unit-cell' : 'cubic',
                       'Wannier_functions_number' : 6,
                       'lattice_parameter' : 4.0184,
                       'highest_occupied_level' : 8.4702,
                       'lowest_unoccupied_level' : 10.5168,
                       'Wannier_file' : 'KTO_hr.dat',
                       'skip_rows' : 150,
                       'manifold' : 't2g',
                       },
              
              'BTB' : {'description' : 'hexagonal BiTeBr, conduction band manifold',
                       'pseudopotential' : 'PAW PBE full-relativistic from PSlibrary',
                       'unit-cell' : 'hexagonal',
                       'Wannier_functions_number' : 6,
                       'lattice_parameter' : 4.2662,
                       'lattice_parameter_c' : 6.487,
                       'highest_occupied_level' : 7.2121,
                       'lowest_unoccupied_level' : 7.5174,
                       'Wannier_file' : 'BTB_hr.dat',
                       'skip_rows' : 576,
                       'manifold' : 'other',
                       },
              
                       
              # 'name' : {'description' : 'Details of the system.',
              #           'pseudopotential' : 'Details of the pseudopotential used in DFT calculations.',
              #           'unit-cell' : 'Unit-cell symmetry (at present, 'cubic' or 'hexagonal').',
              #           'Wannier_functions_number' : 'Number of WFs, it must match the number in the W90 file.',
              #           'lattice_parameter' : 'Lattice parameter a for cubic or hexagonal cell in Angs.',
              #           'lattice_parameter_c' : 'Lattice parameter c in Angs if the unit-cell is hexagonal.',
              #           'highest_occupied_level' : 'Valence band maximum in eV.',
              #           'lowest_unoccupied_level' : 'Conduction band minimum in eV.',
              #           'Wannier_file' : 'Name for the W90 file with extension. Usually name_hr.dat.',
              #           'skip_rows' : 'Number of rows to skip within the W90 file.',
              #           'manifold' : 'Manifold of bands for the system. At present 't2g' or 'other'.',
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

