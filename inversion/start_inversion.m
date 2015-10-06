
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

usr_par.cluster = 'monch';
% 'monch';
% 'euler';
% 'brutus';

usr_par.type = 'source';
% 'source_constrained';
% 'structure';

if( strcmp( usr_par.type, 'structure') )
    [~,~,nx,nz,~,~,~,model_type] = input_parameters();
    [~,rho] = define_material_parameters( nx, nz, model_type);

    usr_par.structure_inversion.rho = rho;
    usr_par.structure_inversion.v0 = 4000;    
end


usr_par.measurement = 'log_amplitude_ratio';
% 'log_amplitude_ratio';
% 'amplitude_difference';
% 'waveform_difference';
% 'cc_time_shift';

% do measurement on displacement or velocity correlations (for NOW: use 'dis')
usr_par.veldis = 'dis';

% filter correlations for inversion; if yes, specify f_min and f_max
usr_par.filter.apply_filter = 'no';
usr_par.filter.f_min = 1/11 - 0.005;
usr_par.filter.f_max = 1/11 + 0.005;

% load array with reference stations and data
usr_par.network = load('../output/interferometry/array_16_ref.mat');
usr_par.data = load('../output/interferometry/data_16_ref_0.mat');

% specify percentile for clipping of kernel      ( = 0 to turn it off )
usr_par.kernel.percentile = 0;

% desing gaussian filter for smoothing of kernel ( = 0 to turn it off )
% usr_par.kernel.smoothing = 0;
usr_par.kernel.smoothing = fspecial('gaussian',[75 75], 30);

usr_par.debug.switch = 'no';
usr_par.debug.df = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inversion parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.tolerance = 1e-3;
options.max_iterations = 1000;
options.init_step_length = 8192.0; 
options.wolfe_try_to_increase_step_length = true;
options.verbose = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% running the inversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);

% start matlabpool and set up path
if( ~strcmp( usr_par.cluster, 'local') )
    addpath(genpath('../'))
    parobj = start_cluster( usr_par.cluster, '', size(usr_par.network.ref_stat,1));
end


% run source inversion
if( strcmp( usr_par.type, 'source' ) )
    
    m0 = make_noise_source();
    m0 = reshape(m0,[],1);
    
    [flag, mfinal, usr_par] = optlib_lbfgs(m0, options, usr_par);

% run source inversion with lower and upper bounds
elseif( strcmp( usr_par.type, 'source_constrained' ) )
    
    type = 'source';
    
    m0 = make_noise_source();
    m0 = reshape(m0,[],1);
    
    xl = 0 * m0;
    xu = inf * m0;
    
    m = projected_steepest_descent(m0,xl,xu,'get_obj_grad',0.05,0);

% run structure inversion
elseif( strcmp( usr_par.type, 'structure' ) )
    
    [~,~,nx,nz] = input_parameters();
    m0 = zeros(nx*nz, 1);
    
    m = steepest_descent(m0,'get_obj_grad',0.05,0);
    m = v0 * ( 1 + m );
    
end


% save solution
save('../output/solution.mat', 'flag', 'mfinal', 'usr_par')


% close matlabpool and clean up path
if( ~strcmp( usr_par.cluster, 'local' ) )
    delete(parobj)
    rmpath(genpath('../'))
end

