
addpath(genpath('../'))
[~,~,nx,nz,~,~,~,model_type,~,n_basis_fct] = input_parameters();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

usr_par.cluster = 'monch';
% 'local';
% 'monch';
% 'euler';
% 'brutus';


usr_par.use_mex = 'yes';


% define v0
[~,rho] = define_material_parameters( nx, nz, model_type);
usr_par.structure_inversion.rho = rho;
usr_par.structure_inversion.v0 = 4000;


% define weighting for kernels
usr_par.kernel.weighting = 0.5;


usr_par.measurement.source = 'log_amplitude_ratio';
usr_par.measurement.structure = 'waveform_difference';
% 'log_amplitude_ratio';
% 'amplitude_difference';
% 'waveform_difference';
% 'cc_time_shift';


% do measurement on displacement or velocity correlations (for NOW: use 'dis')
usr_par.veldis = 'dis';


% filter correlations for inversion; if yes, specify f_min and f_max
usr_par.filter.apply_filter = 'no';
usr_par.filter.f_min = 1/7 - 0.01;
usr_par.filter.f_max = 1/7 + 0.01;


% load array with reference stations and data
usr_par.network = load('../output/interferometry/array_16_ref_small.mat');
usr_par.data = load('../output/interferometry/data_16_ref_0_1h1g_iugg_small.mat');
% usr_par.data = load('../output/interferometry/data_16_ref_0_1h1g_refl1.mat');


% load source distribution that should be used for structure inversion
% (uncomment line if source specified in input_parameters.m should be used)
% usr_par.source_dist = load('initial_models/loga_homog_98.mat');


% specify percentile for clipping of kernel ( = 0 to turn it off )
% be careful: so far not taken into account into gradient
usr_par.kernel.percentile = 0;


% design gaussian filter for smoothing of kernel (set second input variable to [1,1] to turn it off)
usr_par.kernel.imfilter = fspecial('gaussian',[30 30], 12);


% regularization j_total = dj/dm + alpha * ||m-m0||^2  (i.e. set alpha=0 to turn it off)
usr_par.regularization.alpha = 0.007;
usr_par.regularization.weighting = weighting( nx, nz, n_basis_fct );


% debug mode
usr_par.debug.switch = 'no';
usr_par.debug.df = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inversion parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FOR LBFGS
options.tolerance = 1e-3;
options.successive_change = 1e-16;
options.successive_iterations = 100;
options.max_iterations = 100;
options.wolfe_try_to_increase_step_length = false;
options.verbose = true;
options.init_step_length = 32.0;


%%% FOR STEEPEST DESCENT
% options.init_step_length = 1.0;
% options.stepsize = 1.0;
% options.verbose = false;
% options.output_file = 'iterations_lbfgs.tab';
% options.max_iterations = 100;
% options.tolerance = 1e-3;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% running the inversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% start matlabpool and set up path
if( ~strcmp( usr_par.cluster, 'local') )
    addpath(genpath('../'))
    parobj = start_cluster( usr_par.cluster, '', size(usr_par.network.ref_stat,1));
end


% set necessary fields that might not have been set
[usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);


% set up initial model
if( n_basis_fct == 0)
    m_parameters = zeros(nx, nz, 2);
    m_parameters(:,:,1) = make_noise_source();
else
    m_parameters = zeros(nx, nz, n_basis_fct+1);
    m_parameters(:,:,1:n_basis_fct) = make_noise_source();
end
m_parameters(:,:,end) = define_material_parameters(nx,nz,model_type);


% convert to optimization variable
m0 = map_parameters_to_m(m_parameters,usr_par);
usr_par.m0 = m0;


% run inversion
[flag, mfinal, usr_par] = optlib_lbfgs(m0, options, usr_par);


% save solution
save('../output/solution.mat', 'flag', 'mfinal', 'usr_par')


% close matlabpool and clean up path
if( ~strcmp( usr_par.cluster, 'local' ) )
    delete(parobj)
end

