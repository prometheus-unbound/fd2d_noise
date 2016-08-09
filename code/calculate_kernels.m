
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

usr_par.use_mex = 'yes';
% 'yes' (if startup.m indicates a successful compilation)
% 'no' (default)


usr_par.type = 'source';
% 'source'
% 'structure'


usr_par.measurement= 'log_amplitude_ratio';
% 'log_amplitude_ratio';
% 'amplitude_difference';
% 'waveform_difference';
% 'cc_time_shift';


% load data, array and reference stations
usr_par.network = load('../output/array_1_ref.mat');
usr_par.data = load('../output/data_1_ref_model_1_source_gaussian.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate kernels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,~,nx,nz,~,~,~,~,~,~,make_plots] = input_parameters();


% set necessary fields that might not have been set
[usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);


% get source and material
noise_source = make_noise_source('no');
structure = define_material_parameters('no');


% loop over reference stations
n_ref = size(usr_par.network.ref_stat,1);
n_rec = size(usr_par.network.array,1)-1;
t = usr_par.data.t;
nt = length(t);
c_iref = zeros(n_ref,n_rec,nt);

misfit = 0;
gradient = zeros( nx, nz, 2 );

fprintf('\n')

tic
for i_ref = 1:n_ref
        
    
    % each reference station will act as a source once
    src = usr_par.network.ref_stat(i_ref,:);
    rec = usr_par.network.array(  ~ismember( usr_par.network.array, src, 'rows')  , : );
    
    
    % calculate Green function
    if( strcmp( usr_par.type, 'source' ) )
        G_fft = run1_forward_green_mex( structure, src, 0 );
    else
        [G_fft, G] = run1_forward_green_mex( structure, src, 1 );
    end
    
    
    % calculate correlation
    if( strcmp( usr_par.type, 'source' ) )
        c_iref(i_ref,:,:) = run2_forward_correlation_mex( structure, noise_source, G_fft, rec, 0 );
    else
        [c_iref(i_ref,:,:), C] = run2_forward_correlation_mex( structure, noise_source, G_fft, rec, 1 );
    end
    
    
    % calculate misfits and adjoint source functions
    [misfit_iref, adjstf_iref] = make_adjoint_sources( reshape(c_iref(i_ref,:,:),[],nt), usr_par.data.c_data( (i_ref-1)*n_rec + 1 : i_ref*n_rec, : ), t, usr_par.measurement, src, rec );
    
    
    % sum up misfits for all reference stations
    misfit = misfit + sum(misfit_iref);
    
    
    % calculate gradient and sum them up for all reference stations
    if( strcmp( usr_par.type, 'source' ) )
        
        gradient_iref = run3_adjoint( structure, noise_source, G_fft, rec, adjstf_iref, [], 0 );
        gradient = gradient + gradient_iref;
        
    else
        
        [gradient_iref_1, adjoint_state] = run3_adjoint( structure, noise_source, [], rec, adjstf_iref, C, 1 );
        gradient_iref_2 = run3_adjoint( structure, noise_source, [], rec, adjoint_state, G, 0 );
        
        gradient = gradient + gradient_iref_1 + gradient_iref_2;
        
    end
    
    
end
toc

fprintf( 'misfit:   %15.10f\n', misfit )


% reorganize correlation vector
c_iteration = zeros( n_ref*n_rec, length(t) );
for i_ref = 1:n_ref
    c_iteration( (i_ref-1)*n_rec + 1 : i_ref*n_rec, :) = c_iref(i_ref,:,:);
end
toc
