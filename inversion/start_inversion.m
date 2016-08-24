
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

usr_par.type = 'source';
% 'source'
% 'structure'
% 'joint'


usr_par.measurement.type = 'waveform_difference';
% 'log_amplitude_ratio';
% 'amplitude_difference';
% 'waveform_difference';
% 'cc_time_shift';


usr_par.measurement.mode = 'auto';
% manual
% auto


% for design of gaussian filter for smoothing of kernel
usr_par.smoothing.sigma = [1 1];


% provide name of array and of data file
array_file = 'array_nref-1.mat';
data_file = 'correlations_nref-1_model-1_source-gaussian.mat';
% array_file = 'array_nref-1_testing.mat';
% data_file = 'correlations_nref-1_testing.mat';


usr_par.verbose = true;
% true
% false


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inversion parameters - steepest descent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.init_step_length = 1.0;
options.stepsize = 1.0;
options.verbose = true;
options.output_file = 'iterations_steepest_descent.tab';
options.max_iterations = 100;
options.tolerance = 1e-2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% running the inversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- load array and data ----------------------------------------------------
usr_par.network = load([fd2d_path(), 'output', filesep, array_file]);
usr_par.data = load([fd2d_path(), 'output', filesep, data_file]);


%- set necessary fields that might not have been set ----------------------
[usr_par] = usr_par_init_default_parameters(usr_par);


%- set up initial model ---------------------------------------------------
noise_source = make_noise_source('no');
structure = define_material_parameters('no');

m_parameters = zeros(size(structure.mu,1), size(structure.mu,2), 2);
m_parameters(:,:,1) = structure.mu;
m_parameters(:,:,2) = noise_source.distribution;


%- include density and spectrum to usr_par --------------------------------
usr_par.structure.rho = structure.rho;
usr_par.noise_source.spectrum = noise_source.spectrum;


%- convert to optimization variable ---------------------------------------
m0 = map_parameters_to_m(m_parameters, usr_par);


%- run inversion ----------------------------------------------------------
[flag, mfinal, usr_par] = optlib_steepest_descent(m0, options, usr_par);
m_parameters_final = map_m_to_parameters( mfinal, usr_par );


%- plot final model -------------------------------------------------------
plot_models( sqrt( m_parameters_final(:,:,1)  ./ structure.rho), ...
        m_parameters_final(:,:,2), usr_par.network.array, [0, 0, 0, 0]);


%- save solution ----------------------------------------------------------
save([fd2d_path() filesep 'inversion' filesep 'solution.mat'], ...
    'flag', 'mfinal', 'usr_par', '-v7.3')


