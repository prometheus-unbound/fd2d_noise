# fd2d_noise
fd2d_noise is a tool to calculate noise correlations and kernels for various misfit measurements. The forward solution is based on a 2D staggered grid finite difference discretization of the seismic wave equation. The kernels form the basis for an iterative inversion procedure.



INPUT PARAMETERS:
---------------------------------------------------------------------------------------
GENERAL INPUT: /input/input_parameters.m
* size of the computational domain
* number of grid points
* time step and number of time steps
* finite-difference order
* sampling of forward field for adjoint run
* model type as defined in /code/propagation/define_material_parameters.m (be careful, some are hardcoded)
* source type and number of basis functions in frequency domain

ABSORBING BOUNDARY SPECIFICATIONS: /input/absorb_specs.m
* width of boundaries
* on/off switches for each side

OUTPUT SPECIFICATIONS: /input/output_specs.m
* adjoint source path
* verbose switch
* plot switch
* movie switch and location

NOISE SOURCE: /input/make_noise_source.m
* number of noise sources
* spectrum
* spatial distribution

FREQUENCY DOMAIN SPECIFICATIONS: /input/input_interferometry.m
* frequency sampling
* how often the fft for the Green function should be computed



COMPUTING NOISE CORRELATIONS: /code/calculate_data.m
---------------------------------------------------------------------------------------
* specify input in files described above
* specify where the calculation should be done (locally or on Monch/Euler/Brutus)
* define receiver array and which stations will act as reference stations
* results are saved in array_xx_ref.mat and data_xx_ref_yy.mat with xx="number of reference stations" and yy="number of basis functions for the source"


The computation of noise correlations in calculate_data.m basically proceeds in two steps.

STEP 1: Computation of the Green function with source at the reference station
* The source (in earthquake simulations) acts as the reference station in noise correlation simulations.
* The Green function for the reference station is calculated with run_forward1_green.m
* Its Fourier transform is computed on-the-fly and is stored as “G_2_xx.mat” in the directory /output/interferometry/“ with xx="number of reference station"

STEP 2: Computing the actual correlation function
* The calculated Green function together with the specified noise source act as source for the correlation wavefield with run_f



COMPUTING KERNELS: /code/calculate_kernel.m
---------------------------------------------------------------------------------------
* specify input in files described above
* specify kernel type - source or structure
* choose desired measurement for adjoint source calculation
* load respective array_xx_ref.mat and data_xx_ref.mat. for data independent kernels only array is needed; set usr_par.data_independent to 'yes'

The computation of noise source kernels proceeds in three steps.
STEP 1: Computation of "initial" correlations
STEP 2: Measurements and computation of adjoint sources with /tools/make_adjoint_sources.m
STEP 3: Computing the actual kernels with run_noise_source/structure_kernel.m



INVERSION: 
---------------------------------------------------------------------------------------
RUN INVERSION: /inversion/start_inversion.m
* specify input and inversion parameters
* the actual inversion uses scripts by Christian Böhm



CONTRIBUTIONS
---------------------------------------------------------------------------------------
Contributions are always welcome. To contribute to our code, create your own fork, modify it, run git rebase, upload it and send us a pull request.
