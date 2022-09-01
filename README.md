![Logo](https://user-images.githubusercontent.com/94184924/156950182-06b978a4-9302-460f-bb24-5f78c3ecafc8.png)
## BINPO: A CODE FOR ELECTRONIC PROPERTIES OF 2D ELECTRON SYSTEMS

### Version 1.1.0

_BinPo_ is a Python code to compute the electronic band structure and other properties in 2D electron
systems (2DESs). At present, it could be used in cubic and hexagonal systems. For the cubic case, it is possible
to analize the main confinement directions. There is an especial focus in the 2DESs with _t2g_ manifold, like in SrTiO3 (STO) or KTaO3 (KTO),
due to their increasing impact in complex oxide community. So that, in that cases further capabilities are available.
_BinPo_ solves the Schrödinger-Poisson scheme to obtain the self-consistent potential energy along the slab. The 
tight-binding Hamiltonian for this slab is created from the transfer integrals in the Maximally Localized Wannier
Functions (MLWF) basis. Once the self-consistent (SC) solution is found, properties like projected band structure and Fermi
surface, orbital decomposition of electron density and envelope wavefunctions can be computed in post-processing steps.
_BinPo_ gives priority to ease-of-use and efficiency in order to produce realistic simulations at low computational cost.
High quality and customized figures can be obtained from the simulations. To see details of the methodology applied, please
refer to our preprint on https://arxiv.org/abs/2203.11308. For a quick start about how to use this code, please see "Usage" section below.

### Prerequisites

You will need to have installed a version 3.x of Python. If you don't have it, please, refer to the official Python 
website https://www.python.org/. Along with Python 3.x you will need the following libraries:

* NumPy >= 1.19 (https://numpy.org/)
* SciPy >= 1.5 (https://www.scipy.org/)
* Matplotlib >= 3.3 (https://matplotlib.org/)
* NumExpr >= 2.7 (https://numexpr.readthedocs.io/projects/NumExpr3/en/latest/)

We recommend to install a Python scientific distribution like Anaconda (https://www.anaconda.com/) or another one, in 
which the modules above mentioned are already included. Additionally, you will have to install the Atomic Simulation
Enviroment (ASE) library, version >= 3.2 (https://wiki.fysik.dtu.dk/ase/).

### Features

Currently, _BinPo_ can be used for cubic and hexagonal systems. For the cubic case, a rotation algorithm allows you to
analyze the main confinement directions ((001), (110) and (111)). On the other hand, for hexagonal systems the confinement
direction of the 2DES must be along the c axis. We already included several W90 files in _BinPo_ for the following systems:
* strontium titanate, STO
* potassium tantalate, KTO
* BiTeBr

Additionally, the code has modularity to append new materials (see _~/BPexamples/Adding_W90files.pdf_).

_BinPo_ allows for computing the following properties:
* The SC solution to the Schrödinger-Poisson scheme in the slab. It includes the SC-potential energy, electron density,
  electric field and relative permittivity values. 
* Bandstructure of the 2DES along high-symmetry points within the irreducible Brillouin zone. The bands calculation
  can be total (without projections) or projected onto a set of planes. In the particular case of _t2g_ 2DESs, projections
  onto the atomic orbitals are also available.
* Fermi surface and, more generally, energy slices. In the particular case of _t2g_ 2DESs, projections onto the atomic
  orbitals are also available.
* For _t2g_ 2DESs there are two more additional capabilities: decomposition of electron density into the atomic
  orbitals contributions and calculation of envelope wavefunctions at gamma point. 

Furthermore, you can analyze the simulations results by Matplotlib interactive plots and generate customized high
quality plots for publishing.

### License

_BinPo_ is Copyrighted by (C) 2021 BinPo Team. This code is distributed under the terms of the GNU General Public 
License v3 (GPLv3), see ~/COPYING file or http://www.gnu.org/copyleft/gpl.txt. So that, you are free to use, modify,
improve and redistribute it.

### Contact and help

For reporting bugs or asking questions about the code, please write to the following email: emanuelm@ucm.es.
Request for adding features will be welcomed but without guarantees of future implementations.

### What does each file/folder mean?

* _WFolder_:         This folder holds the initial Wannier90 output files (W90 files). These files come from a 
                   DFT calculation + Wannierization, and can be seen as the real space 
                   Hamiltonian of the bulk unit cell. There is a specific file for each material, bands manifold,
                   as well as for the characteristics of the original DFT calculation (e. g. XC-functional).

* _BPexamples_:      This folder contains _.pdf_ files with detailed descriptions about the calculation presented as 
                   examples in our main manuscript. Also a stepwise guide to add new W90 files can be found.

* _config_files_:    This folder contains the configuration files _.yaml_ for the _BinPo_ components. Except for _BP-preproc.py_
		   and _BP-fast_plot.py_ components (which not require a configuration file), each of the _BinPo_ components, named
		   as _BP-components.py_, has an associated _component.yaml_ configuration file, where by default the values not set 
		   by terminal are taken from. Besides, in this folder there is a _help_config.md_ file which explains the meaning of
		   all parameters in the _.yaml_ files.

* _BP-preproc.py_:   Pre-processing component. It performs the tasks of filtering and rearrangement from the initial
                   Wannier files.

* _BPmodule.py_:     This file is the main module and contains classes, methods an attributes.
                   It is only called by others programs, the file itself does not need to be executed.	

* _BPdatabase.py_:   This module contains a dictionary with information about the materials whose W90 files are included.
                   More data could be added if the corresponding W90 files are available. This file is not runnable, because
                   it is imported by other programs and the information is obtained through the 'material' keyword.
		
* _BP-scp.py_:       Component for calculation of the SC potential energy along the slab (SC-potential).
 				
* _BP-fast_plot.py_: Quick plotter for the output of _BP-scp.py_.

* _BP-bands.py_:     Band structure calculator component. Bands calculation can include projections onto planes. For the 
		   particular case of _t2g_ 2DESs, projections onto the orbital character is also available.
                   You can pass the identifier for band calculation by terminal, and many other parameters.
 
* _BP-energy_slices.py_:  
		   Energy slices calculator component. The calculation can include orbital character for _t2g_ 2DESs.
		   You can pass the identifier for slice calculation by terminal, and many other parameters. At the end,
		   the file is automatically saved for plotting.

* _BP-energy_plot.py_:
      		   Energy slices plotter component. Plot the output of _BP-energy_slices.py_.

* _BP-envelope_wfs.py_: 
	           At present, it is just available for _t2g_ 2DESs. Envelope wavefunctions calculator at gamma point.
	           It allows for obtaining and visualizing the envelope wavefunctions around the gamma point.

* _BP-orb_density.py_: 
	           At present, it is just available for _t2g_ 2DESs. Orbital decomposition of electron density.
	           It allows for decomposing and plotting the electron density according to the orbital character.
 
* _COPYING.md_:      GPLv3 license file.

* _README.md_:       This file.


### Usage

To see the usage and updatable parameters from command-line, please type:

     $ python BP-component.py -h

where _BP-component.py_ is any _BinPo_ component. A list of the more frequently adjustable parameters will appear.
If a parameter is ommited, its value will be taken by default from the corresponding configuration file.

In the following lines you will find a short description of how to use _BinPo_. For further details you can take a
look at _~/BPexamples_. 

###    STEP 1: pre-processing step:

Run "BP-preproc.py". This component has two mandatory arguments: material (mt) and confinement direction (cfd). Choose the combination you want,
for example, STO and 111, as follows:

	$ python BP-preproc.py -mt STO -cfd 111

It will generate a folder, named as "Hr" + "material" + "direction", that holds the files corresponding to the rearrangment and
filtering of r-space Hamiltonian, discretized along normal direction according to the number of planes. You will find the folder
 _~/HrSTO111_ after a succesful pre-processing.
			
**NOTE1:** It should be noticed that you need to execute this step just once for each W90 file/material/direction combination.

**NOTE2:** For hexagonal systems rotations are not allowed, so that the (001) is the confinement direction. 


###    STEP 2: SC-potential calculation:
		
Run _BP-scp.py_ component with a defined job identifier (_id_) as: 

	$ python BP-scp.py -id identifier -mt STO -cfd 111
	
In this simplified example, many other parameters are being ignored. These will be taken from _~/conf_files/scp.yaml file_. 
The identifier will be unique for this calculation and can be recall later in any post-processing steps. After a successfull SC calculation,
you will obtain an output folder holding the SC solution, a _.log_ file and a _.yaml_ file with a dictionary of parameters used in the 
calculation. You can quickly check the output file by means of _BP-fast_plot.py_ as:
     
	 $ python BP-fast_plot.py -id identifier
				
###    STEP 3: post-processing step:

**NOTE:** Each of the post-processing components needs a previous SC-potential calculation whose solution will be accessed by the specific 
identifier.

Recall the SC solutions and their features from any of the post-processing components (_BP-bands.py_, 
_BP-energy_slices.py_, etc.). Check the parameters needed in each particular case. All new runs will be automatically logged in the 
file _~/identifier/identifier.log_. For example, to compute the bandstructure along KGK path (_ph_), run the line below. In this case, 
the ommited parameters will be taken from _~/conf_files/bands.yaml_ file.

	$ python BP-bands.py -id identifier -ph KGK

### Further information:

* For learning about the different options to edit plots (colors, colormaps, formats, etc.), check the Matplotlib documentation
  at https://matplotlib.org/.
* For general information about .yaml files, visit https://yaml.org/.

