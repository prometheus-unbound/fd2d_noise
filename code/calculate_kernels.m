
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

usr_par.use_mex = 'no';
% 'yes' (if startup.m indicates a successful compilation)
% 'no' (default)


usr_par.type = 'structure';
% 'source'
% 'structure'


usr_par.measurement= 'waveform_difference';
% 'log_amplitude_ratio';
% 'amplitude_difference';
% 'waveform_difference';
% 'cc_time_shift';


% provide name of array and of data file
array_file = 'array_1_ref.mat';
data_file = 'correlations_1_ref_model_1_source_gaussian.mat';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate kernels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check path and mex files
fd2d_path = fd2d_path();
check_mex_files( usr_par.use_mex );  


[~,~,nx,nz,dt,nt,~,model_type,source_type,~,make_plots] = input_parameters();


% load data 
usr_par.network = load([fd2d_path() 'output' filesep array_file]);
usr_par.data = load([fd2d_path() 'output' filesep data_file]);


% set necessary fields that might not have been set
[usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);


% get source and material
noise_source = make_noise_source('no');
structure = define_material_parameters('no');


% loop over reference stations
n_ref = size(usr_par.network.ref_stat,1);
n_rec = size(usr_par.network.array,1)-1;
t = -(nt-1)*dt:dt:(nt-1)*dt;
nt = length(t);

if( ~exist( filename('correlations', n_ref), 'file') )
    correlations = zeros(n_ref,n_rec,nt);
else
    correlations = parload( filename('correlations', n_ref) );
end

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
        
        if( ~exist( filename('G_fft', i_ref), 'file' ) )
            G_fft = run1_forward_green_mex( structure, src, 0 );
            parsave( filename('G_fft', i_ref), G_fft, [] )
        else
            G_fft = parload( filename('G_fft', i_ref) );
        end
        
    else
        [G_fft, G] = run1_forward_green_mex( structure, src, 1 );
    end
    
    
    % calculate correlation
    if( strcmp( usr_par.type, 'source' ) )
        
        if( ~exist( filename('correlations', n_ref), 'file') )
            correlations(i_ref,:,:) = run2_forward_correlation_mex( structure, noise_source, G_fft, src, rec, 0 );
        end
        
    else
        [correlations(i_ref,:,:), C] = run2_forward_correlation_mex( structure, noise_source, G_fft, src, rec, 1 );
    end
    
    
    % calculate misfits and adjoint source functions
    [misfit_iref, adjstf_iref] = make_adjoint_sources( reshape(correlations(i_ref,:,:),[],nt), reshape(usr_par.data.correlations( i_ref,:,:),[],nt), t, usr_par.measurement, src, rec );
    
    
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
if( ~exist( filename('correlations', n_ref), 'file') )
    save( filename('correlations', n_ref), correlations, [])
end


%% plot kernel missing
if( strcmp( usr_par.type, 'source' ) )
    plot_models( [], gradient(:,:,2), usr_par.network.array, [0 0 0 0]);
else
    plot_models( gradient(:,:,1), [], usr_par.network.array, [0 0 0 0]);
end

