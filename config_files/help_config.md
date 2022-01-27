MEANING OF THE PARAMETERS IN THE CONFIGURATION FILES
====================================================


In the following lines you will find the specifications for the .yaml configuration files in BinPo.

* Please, refer to the official Python website to check syntaxis if needed. 
* For learning about the differents option to edit plots (colors, colormaps, formats, etc.),
check the Matplotlib documentation at https://matplotlib.org/.
* For general information about .yaml files, visit https://yaml.org/.


-----------------------------------------------------------------------------------------------------------------
 scp.yaml
-----------------------------------------------------------------------------------------------------------------

SCP_CALCULATION:
---------------

- identifier = String. Identifier for the calculation. It is unique and can be recalled later in post-processing.
- material = String. Name of the material (cubic perovskite ABO3) to be used. Allowed names: STO, KTO.
- crystal_face = String. Crystal face expressed as string 'hkl'. Allowed values: '100' ('010','001'), '110' 
  ('011','101') and '111'.
- number_of_planes = Integer. Total number of planes stacked along normal direction. Allowed range: [10, 200].
- shift_from_LUL = Float. Energy shift (dE) from LUL in eV. Used to define the Fermi level as ef = LUL + dE. 
  Allowed range: [0.0, 0.02].
- BC1_topmost_layer = Float. Dirichlet boundary condition for V in eV at the topmost layer (i.e. V[0]). Allowed 
  range: [-1, -0.001]. 
- BC2_in_bulk = Float. Dirichlet boundary condition for V in eV in the bulk (i.e. V[L-1]). Allowed range: 
  [bc1, 3*|bc1|].
- Neumann_at_bulk = Boolean. Whether or not applied Neumann boundary condition at V[L-1]. Overrides the Dirichlet 
  option for BC2_in_bulk.
- sqrt_kgrid_numbers = Integer. Square root of the total number of kpoints in k-grid. The total 2D Brillouin zone 
  is spanned. It must be > 10.
- k_shift = Float (list). 2D kgrid offset as follows: (k_x, k_y) + (k_shift[0], k_shift[1]).
- temperature = Float. Temperature in K. Allowed range: [0.001, 315.0]
- mixing_factor = Float. Mixing factor for the over-relaxation mixing method. Tipically 0.06-0.4.
- permittivity_model = String. Model used for relative permittivity as function of electric field E. Pre-defined 
  models are: 'Cop', 'Ang' and 'cte_N'. Alternatively, an expression as function of E in python syntax could be 
  introduced. For a quick explanation of the models see below.* 
- potential_live_visualization = Boolean. Step by step visualization of the Vin and Vout in the iterative process.
- error_live_visualization = Boolean. Step by step visualization of the total error as funcion of iterations.

- ADVANCED_PARAMETERS:
  -------------------
  
  - conv_threshold = Float. Convergence threshold for V in self-consistent calculation.
  - max_iterations = Integer. Maximum number of allowed iterations.
  - Total_Hk_method = String. Method by which total Hamiltonian is constructed. It can be either 'vectorized' or 
    'iterable'. 'vectorized' option is suggested whenever it is possible.
  - V_initial = String. Kind of z dependency of initial potential energy. It can be either 'linear', 
    V(z) = z*(BC2-BC1)/(L-1) + BC1 or 'exponential', V(z)  = BC1*exp(-z/2), where BC1, BC2 are the boundary conditions
    and L is the number_of_planes.
  - cons_charge_background = Boolean. Whether or not to include a constant charge in Poisson equation.
  - charge_per_site = Float. If 'cons_charge_background = yes', it sets the value of the charge per site.
  - charge_extension = Integer. If 'cons_charge_background = yes', it sets the extension of charge_per_site from 
    top-most layer to the bulk. It must be < (number_of_planes-1).


-----------------------------------------------------------------------------------------------------------------
 bands.yaml
-----------------------------------------------------------------------------------------------------------------

BAND_STRUCTURE:
---------------

- identifier = String. Identifier for the calculation.
- path = String. Path in the BZ1 to perform the bandstructure calculation. Special points for (100) direction are G,X 
  and M, for (110) direction are G,X,Y,S and for (111) direction are G,K,M.
- number_of_kpoints = Integer. Number of k-points along the path.
- reference_kpoint = String. Special point used to set the 0 in k axis. Tipically "G" point.             
- Total_Hk_method = String. Method by which total Hamiltonian is constructed. It can be either 'vectorized' or 'iterable'.
  'vectorized' option is suggested whenever it is possible.
- num_bands = Integer. Number of bands to compute. It must be between 1 and 6*number_of_planes.
- bands_task = Integer. Task to perform. Options are 0:'simple_bandstructure', 1:'orbital_character', 2:'plane_projection'
  and 3:'all'.
- initial_plane = Integer. Initial plane from which the projection will start. Only used if bands_task == 2 or == 3.
- final_plane = Integer. Upper limit for the projections. Planes involved will be (plane_init, plane_fin-1). Only used if 
  bands_task == 2 or == 3.

- TOTAL_BANDS:
--------------
             
  - PLOT_ADJUST:
----------------
   - plotstyle = String. Matplotlib plot style. 
   - xy_limits = Float (list). Limits for the plot. Write it as [x_min, x_max, y_min, y_max].
   - linecolor = String. Matplotlib color for the lineplot.                             
   - linewidth = Float. Linewidth.
   - fig_size = Float (list). Size of the matplotlib figure set as [x_dim, y_dim].
   - axis_adjust = Float (list). Axes location within the figure set as [left, bottom, right, top].
   - title = String. Title for the figure. If title == "", it will be ignored.
   - title_size = Float. Fontsize for the figure title. Only if title != "".
   - shadow_above_Ef = Float. Opacity above Fermi level. It must be between 0.0 (no opacity) and 1.0 (solid).
   - shadow_color = String. Matplotlib color for the opacity above Fermi level.

  - LABELS: 
-----------   
   - xlabel = String. X label for the figure. It accepts LaTeX syntax.
   - xfontsize = Float. Fontsize for the x label.
   - ylabel = String. Y label for the figure. It accepts LaTeX syntax.
   - yfontsize = Float. Fontsize for the y label.
   - ticksize = Float. Size of the ticks.

  - SAVING: 
-----------
   - save_bands = Boolean. Whether or not to save the bands to a file.
   - save_plot = Boolean. Whether or not to save the plot to a file.
   - format = String. Format for the plot if saved.
   - dpi = Integer. Dots per inch for the plot if saved. 

- ORBITAL_CHARACTER:
--------------------

  - PLOT_ADJUST:
----------------
   - plotstyle = String. Matplotlib plot style. 
   - xy_limits = Float (list). Limits for the plot. Write it as [x_min, x_max, y_min, y_max].
   - color_seq = String. Sequence of three colors as color1,color2,color3.                            
   - point_size = Float. Pointsize for the scatter plot.
   - fig_size = Float (list). Size of the matplotlib figure set as [x_dim, y_dim].
   - axis_adjust = Float (list). Axes location within the figure set as [left, bottom, right, top].
   - title = String. Title for the figure. If title == "", it will be ignored.
   - title_size = Float. Fontsize for the figure title. Only if title != "".
   - shadow_above_Ef = Float. Opacity above Fermi level. It must be between 0.0 (no opacity) and 1.0 (solid).
   - shadow_color = String. Matplotlib color for the opacity above Fermi level.
             
  - COLOR_TRIANGLE:
-------------------
   - proportion = String. Axes percentage occupied by the color triangle.
   - location = Integer. Location according to matplotlib positions rules. Values are between [1,10]
   - padding = Float. Separation from the axes if triangle is near to some edge.
   - fontsize = Float. Fontsize for text in triangle color.
    
  - LABELS:
-----------                   
   - xlabel = String. X label for the figure. It accepts LaTex syntax.
   - xfontsize = Float. Fontsize for the x label.
   - ylabel = String. Y label for the figure. It accepts LaTeX syntax.
   - yfontsize = Float. Fontsize for the y label.
   - ticksize = Float. Size of the ticks.
     
  - SAVING: 
-----------
   - save_bands = Boolean. Whether or not to save the bands to a file.
   - save_plot = Boolean. Whether or not to save the plot to a file.
   - format = String. Format for the plot if saved.
   - dpi = Integer. Dots per inch for the plot if saved. 

- PLANE_PROJECTION: 
-------------------
  
  - PLOT_ADJUST:
  --------------
   - plotstyle = String. Matplotlib plot style.  
   - xy_limits = Float (list). Limits for the plot. Write it as [x_min, x_max, y_min, y_max].
   - colormap = String. Matplotlib colormap.
   - background_color = String. Matplotlib color for the background.
   - point_size = Float. Pointsize for the scatter plot. 
   - fig_size = Float (list). Size of the matplotlib figure set as [x_dim, y_dim].
   - axis_adjust = Float (list). Axes location within the figure set as [left, bottom, right, top].
   - title = String. Title for the figure. If title == "", it will be ignored.
   - title_size = Float. Fontsize for the figure title. Only if title != "".
   - shadow_above_Ef = Float. Opacity above Fermi level. It must be between 0.0 (no opacity) and 1.0 (solid).
   - shadow_color = String. Matplotlib color for the opacity above Fermi level.
     
  - COLORBAR:
-------------  
   - location = Float (list). Location of the colorbar set as [x, y, width, height]
   - textbar = String (list). Text to be located at the bottom and at the top of the colorbar.
   - fontsize  = Float. Fontsize for the text in colorbar.
   - fontcolor = String. Color for the text in colorbar.
    
  - LABELS:
----------    
   - xlabel = String. X label for the figure. It accepts LaTex syntax.
   - xfontsize = Float. Fontsize for the x label.
   - ylabel = String. Y label for the figure. It accepts LaTex syntax.
   - yfontsize = Float. Fontsize for the y label.
   - ticksize = Float. Size of the ticks.
   
  - SAVING: 
-----------
   - save_bands = Boolean. Whether or not to save the bands to a file.
   - save_plot = Boolean. Whether or not to save the plot to a file.
   - format = String. Format for the plot if saved.
   - dpi = Integer. Dots per inch for the plot if saved. 
    
 
-----------------------------------------------------------------------------------------------------------------
 energy_slices.yaml
----------------------------------------------------------------------------------------------------------------- 
 
ENERGY_SLICES:
--------------
           
- identifier = String. Identifier to compute the energy slice.
- sqrt_kgrid_numbers = Integer. Square root of the total number of points in k-grid. The total 2D Brillouin zone
  is spanned.
- kbox_factor = Float. Factor to modify the box limits in k-space. Value 1.0 corresponds to the BZ1.
- kbox_shift = Float (list). 2D k-grid offset as follows: (k_x, k_y)*kbox_factor + (kbox_shift[0], kbox_shift[1]).
- win_energy_calc = Float. Energy window to get the eigenvalues in SciPy routine of diagonalization.
- batches = Integer. Number of batches to split the total k-points. K-points per batch will be 'sqrt_kgrid_numbers**2/batches'.
  In consequence, 'batches' has to be divisor of 'sqrt_kgrid_numbers**2', otherwise an error will arise.
- energy_cut = Float. Energy in eV at which the slice is computed. 0.0 is the Fermi level.
- outfile = String. Name of the output file. If 'default' it will be saved as 'identifier_ES.dat'

-----------------------------------------------------------------------------------------------------------------
 energy_plot.yaml
----------------------------------------------------------------------------------------------------------------- 

ENERGY_PLOTTER:
---------------
           
- identifier = String. Identifier to plot the energy slice.
- input_file = String. Name of the input file. If 'default' it will load the deafault name 'identifier_ES.dat'.
- orbital_char = Boolean. Whether or not to plot including orbital character.
- color_seq = String. Sequence of 3 colors as color1,color2,color3.  
- energy_window = Float. Energy window. The points shown will range from (energy_cut-energy_window) to energy_cut.
                   
- PLOT_ADJUST:
-------------- 
 - plotstyle = String. Matplotlib plot style. 
 - point_size = Float. Pointsize for the scatter plot.
 - xy_limits = Float (list). Limits for the plot. Write it as [x_min, x_max, y_min, y_max].
 - fig_size = Float (list). Size of the matplotlib figure set as [x_dim, y_dim].
 - axis_adjust = Float (list). Axes location within the figure set as [left, bottom, right, top].
 - title = String. Title for the figure. If title == "", it will be ignored.
 - title_size = Float. Fontsize for the figure title. Only if title != "".
   
- LABELS:
---------
  - xlabel = String. X label for the figure. It accepts LaTex syntax.
  - xfontsize = Float. Fontsize for the x label.
  - ylabel = String. Y label for the figure. It accepts LaTex syntax.
  - yfontsize = Float. Fontsize for the y label.
  - ticksize = Float. Size of the ticks.
   
- COLOR_TRIANGLE:
-----------------
  - proportion = String. Axes percentage occupied by the color triangle.
  - location = Integer. Location according to matplotlib positions rules. Values are between [1,10]
  - padding = Float. Separation from the axes if triangle is near to some edge.
  - fontsize = Float. Fontsize for text in triangle color.
   
- SAVING:
---------
  - save_plot = Boolean. Whether or not to save the plot to a file.
  - format = String. Format for the plot if saved.
  - dpi = Integer. Dots per inch for the plot if saved.  

 
-----------------------------------------------------------------------------------------------------------------
 envelope_wfs.yaml
----------------------------------------------------------------------------------------------------------------- 
 
ENVELOPE_WAVEFUNCTIONS:
-----------------------

- identifier = String. Identifier to compute the envelope wavefunctions at gamma point.
- intensity_factor = Float. Scale factor to regulate the intensity of the envelope wavefunctions.                  
- number_of_wavefunctions = Integer. Number of envelope wavefunctions to be computed.
          
- PLOT_ADJUST:
--------------
  - plotstyle = String. Matplotlib plot style. 
  - linewidth_wfs = String. Linewidth for the wavefunctions curves. 
  - linewidth_V = String. Linewidth for the potential profile. 
  - markersize = Float. Scatter size for the potential profile.
  - axis_adjust = Float (list). Axes location within the figure set as [left, bottom, right, top].
  - fig_size = Float (list). Size of the matplotlib figure set as [x_dim, y_dim].
  - color_V = String. Color for the potential profile.
  - color_wfs = String. Color for the envelope wfs curves.
  - alpha = Float. Opacity of the envelope wfs. curves. It must be between 0.0 (transparent) and 1.0 (solid).
  - title = String. Title for the figure. If title == "", it will be ignored.
  - title_size = Float. Fontsize for the figure title. Only if title != "".
     
- LABELS:
---------
  - xlabel = String. X label for the figure. It accepts LaTex syntax.
  - xfontsize = Float. Fontsize for the x label.
  - ylabel = String. Y label for the figure. It accepts LaTex syntax.
  - yfontsize = Float. Fontsize for the y label.
  - ticksize = Float. Size of the ticks.
  - legends = Boolean. Whether or not to include the legends in the plot.
  - legend_size = Float. If 'legends = yes', it sets the fontsize of legends.
   
- SAVING:
---------
  - save_data = Boolean. Whether or not to save the data to a file.
  - save_plot = Boolean. Whether or not to save the plot to a file.
  - format = String. Format for the plot if saved.
  - dpi = Integer. Dots per inch for the plot if saved.  
 
 
-----------------------------------------------------------------------------------------------------------------
 orb_density.yaml
-----------------------------------------------------------------------------------------------------------------  
 
DENSITY_ANALYSIS:
-----------------

- identifier = String. Identifier to compute the orbital decomposition of the electron density.
- sqrt_kgrid_numbers = Square root of the total number of points in the k-grid. The total 2D Brillouin zone
  is spanned. It must be > 10.
- k_shift = Float (list). 2D kgrid offset as follows: (k_x, k_y) + (k_shift[0], k_shift[1]).
     
- PLOT_ADJUST:
--------------
  - plotstyle = String. Matplotlib plot style.
  - linewidth = Float. Linewidth for all profiles.
  - markersize = Float. Scatter size for all profiles.
  - axis_adjust = Float (list). Axes location within the figure set as [left, bottom, right, top].
  - fig_size = Float (list). Size of the matplotlib figure set as [x_dim, y_dim].
  - color_seq = String. Sequence of 4 colors as color1,color2,color3,color4.
  - fill_curves = Boolean. Whether or not to fill the area under the curves.
  - alpha = Float. If 'fill_curves = yes', it sets the opacity of the encompassed area. 
  - title = String. Title for the figure. If title == "", it will be ignored.
  - title_size = Float. Fontsize for the figure title. Only if title != "".
   
- LABELS:
---------
  - xlabel = String. X label for the figure. It accepts LaTex syntax.
  - xfontsize = Float. Fontsize for the x label.
  - ylabel = String. Y label for the figure. It accepts LaTex syntax.
  - yfontsize = Float. Fontsize for the y label.
  - ticksize = Float. Size of the ticks.
  - legends = Boolean. Whether or not to include the legends in the plot.
  - legend_size = Float. If 'legends = yes', it sets the fontsize of legends.
    
- SAVING:
---------
  - save_data = Boolean. Whether or not to save the data to a file.
  - save_plot = Boolean. Whether or not to save the plot to a file.
  - format = String. Format for the plot if saved.
  - dpi = Integer. Dots per inch for the plot if saved.  
 
 
