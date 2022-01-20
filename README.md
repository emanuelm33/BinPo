## BINPO: A CODE FOR ELECTRONIC PROPERTIES OF 2D ELECTRON SYSTEMS

### Version 1.0
<div style="text-align: justify">
BinPo is a Python code to compute the electronic bandstructure and other properties in 2D electron
systems (2DES). At present, it could be used for electronic properties at the surface/interface of perovskite  
complex oxides, like strontium titanate (STO), in the main crystal faces. BinPo solves the Schrödinger-Poisson 
scheme to obtain the self-consistent potential energy along the slab. The tight-binding Hamiltonian for this slab 
is created from the transfer integrals in the Maximally Localized Wannier Functions (MLWF) basis. Once the self-
consistent (SC) solution is found, properties like projected bandstructure and Fermi surface, orbital decomposition
of electron density and envelope wavefunctions can be computed in post-processing steps. BinPo gives priority to 
ease-of-use and efficiency in order to produce realistic simulations at low computational cost. High quality and 
customized figures can be obtained from the simulations. To see details of the methodology applied, please refer to 
our manuscript. To know how to use this code, please see "Usage" section below.
</div>
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

Currently, BinPo can be used for 2DES formed on the (100), (110) and (111) faces of:
* strontium titanate, STO
* potassium tantalate, KTO
Additionally, the code has modularity to append new materials (cubic ABO3 perovskites) if desired. 
BinPo allows for computing the following properties:
* The SC solution to the Schrödinger-Poisson scheme in the slab. It includes the SC-potential energy, electron density,
  electric field and relative permittivity values. 
* Bandstructure of the 2DES along high-symmetry points within the irreducible Brillouin zone. The bands calculation
  can be total (without projections), projected onto the atomic orbitals or projected onto a set of planes.
* Fermi surface and, more generally, energy slices with or without projections onto the atomic orbitals.
* Decomposition of electron density into the atomic orbitals contributions.
* Envelope wavefunctions at gamma point calculation. 
Furthermore, you can analyze the simulations results by Matplotlib interactive plots and generate customized high
quality plots for publishing.

### License

BinPo is Copyrighted by (C) 2021 BinPo Team. This code is distributed under the terms of the GNU General Public 
License v3 (GPLv3), see ~/COPYING file or http://www.gnu.org/copyleft/gpl.txt. So that, you are free to use, modify,
improve and redistribute it.

### Contact and help

For reporting bugs or asking questions about the code, please write to the following email: emanuelm@ucm.es.
Request for adding features will be welcomed but without guarantees of future implementations.

### What does each file/folder mean?

* WFolder:         This folder holds the initial Wannier90 output files (W90 files). These files come from a 
                   DFT calculation + Wannier interpolation, and can be seen as the real space 
                   Hamiltonian of the bulk unit cell. There is a specific file for each material,
                   as well as for the characteristics of the original DFT calculation (e. g. XC-functional).

* examples:        This folder contains .pdf files with detailed descriptions about the calculation presented as 
                   examples in our main manuscript.

* config_files:    This folder contains the configuration files .yaml for each component of BinPo. Besides, there
                   is a "help_config.txt" file which explains the meaning of all parameters.

* BP-preproc.py:   Preprocessing component. It performs the tasks of filtering and rearrangement from the initial
                   Wannier files.

* BPmodule.py:     This file is the main module and contains classes, methods an attributes.
                   It is only called by others programs, the file itself does not need to be executed.	

* BPdatabase.py:   This module contains a dictionary with the different materials and DFT treatments info.
                   For the time being, it contains the info for STO and KTO calculations. More data could be
                   added if the corresponding Wannier files are available. This file is not runnable, because
                   is imported by other programs and the information is obtained through the 'material' keyword.
		
* BP-scp.py:       Component for calculation of the SC potential energy along the slab (SC-potential).
                   Parameters not passed through the terminal are set by default from ~/conf_files/scp.yaml.
				
* BP-fast_plot.py: Quick plotter for the output of "BP-scp.py".

* BP-bands.py:     Bandstructure calculator component. Bands calculation can include orbital character and/or bands 
		   projections onto planes. You can pass the identifier for band calculation by terminal, and many 
		   other parameters. Parameters not set will be taken from ~/conf_files/bands.yaml

* BP-energy_slices.py:  
		   Energy slices calculator component. The calculation includes orbital character. You can pass the 
		   identifier for slice calculation by terminal, and many other parameters. Parameters not set will
	           be taken form ~/conf_files/energy_slices.yaml. At the end, the file is automatically saved for 
	           plotting.

* BP-energy_plot.py:
      		   Energy slices plotter component. Plot the output of BP-energy_slices.py. Parameters not passed through 
	           the terminal will be set by default from ~/conf_files/energy_plot.yaml.

* BP-envelope_wfs.py: 
	           Envelope wavefunctions calculator at gamma point. It allows for obtaining and visualizing the envelope 
	           wavefunctions around the gamma point. Parameters not passed through the terminal will be set by default
	           from ~/conf_files/envelope_wfs.yaml.

* BP-orb_density.py: 
	           Orbital decomposition of electron density. It allows for decomposing and plotting the electron density 
	           according to the orbital character. Parameters not passed through the terminal will be set by default
	           from ~/conf_files/orb_density.yaml.

* COPYING.md:      GPLv3 license file.

* README.md:       This file.


### Usage

To see the usage and updatable parameters from command-line, please type:
     $ python BP-component.py -h
where BP-component.py is any BinPo component. A list of the more frequently adjustable parameters will appear.
If a parameter is ommited, its value will be taken by default from the corresponding configuration file.

In the following lines you will find a short description of how to use BinPo. For further details you can take a
look at ~/BPexamples folder. 

###    STEP 1: pre-processing step:

Run "BP-preproc.py". This component has two mandatory arguments: material (mt) and confinement direction (cfd). Choose the combination you want,
for example, STO and 111, as follows:

	$ python BP-preproc.py -mt STO -cfd 111

It will generate a folder, named as "Hr" + "material" + "direction", that holds the files corresponding to the rearrangment and
filtering of r-space Hamiltonian, discretized along normal direction according to the number of planes. You will find the folder
 ~/HrSTO111 after a succesful pre-processing.
			
NOTE: It should be noticed that you need to execute this step just once for each W90 file/material/direction combination.

###    STEP 2: SC-potential calculation:
		
Run "BP-scp.py" component with a defined job identifier (id) as: 

	$ python BP-scp.py -id identifier
	
In this simplified example, many other parameters are being ignored. These will be taken from ~/conf_files/scp.yaml file. 
The identifier will be unique for this calculation and can be recall later in any post-processing step. If the calculation converges
you will obtain an output folder holding the SC solution, a .log file and a .yaml file with a dictionary of parameters used in the 
calculation. You can quickly check the output file by means of "BP-fast_plot.py" as:
     
	 $ python BP-fast_plot.py -id identifier
				
###    STEP 3: post-processing step:

NOTE: Each of the post-processing components needs a previous SC-potential calculation whose solution will be accessed by the specific 
identifier.

Recall the SC solutions and their features from any of the post-processing components (BP-bands.py, 
BP-energy_slices.py, etc.). Check the parameters needed in each particular case. All new runs will be automatically logged in the 
file ~/identifier/identifier.log. For example, to compute the bandstructure along KGK path (ph), run the line below. In this case, 
the ommited parameters will be taken from ~/conf_files/bands.yaml file.

	$ python BP-bands.py -id identifier -ph KGK

### Further information:

* For learning about the differents option to edit plots (colors, colormaps, formats, etc.), check the Matplotlib documentation
  at https://matplotlib.org/.
* For general information about .yaml files, visit https://yaml.org/.

