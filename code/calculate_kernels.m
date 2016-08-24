
%==========================================================================
% user input
%==========================================================================

usr_par.type = 'source';
% 'source'
% 'structure'


usr_par.measurement.type = 'log_amplitude_ratio';
% 'log_amplitude_ratio';
% 'amplitude_difference';
% 'waveform_difference';
% 'cc_time_shift';


usr_par.measurement.mode = 'manual';
% manual
% auto


usr_par.data_independent = 'no';
% yes
% no


% provide name of array and of data file
array_file = 'array_nref-1.mat';
data_file = 'correlations_nref-1_model-1_source-gaussian.mat';


usr_par.verbose = false;
% true
% false


%==========================================================================
% calculate kernels
%==========================================================================

%- check path -------------------------------------------------------------
fd2d_path = fd2d_path();


%- get configuration and set up time vector -------------------------------
[~, ~, nx, nz, dt, nt, ~, model_type, source_type, ~, make_plots] = input_parameters();
t = - (nt - 1) * dt:dt:(nt - 1) * dt;
nt = length(t);


%- load array -------------------------------------------------------------
usr_par.network = load([fd2d_path(), 'output', filesep, array_file]);
n_ref = size(usr_par.network.ref_stat, 1);
n_rec = size(usr_par.network.array, 1) - 1;


%- load data --------------------------------------------------------------
if (strcmp(usr_par.data_independent, 'no'))
    usr_par.data = load([fd2d_path(), 'output', filesep, data_file]);
else
    usr_par.data.correlations = zeros(n_ref, n_rec, nt);
end


%- set necessary fields that might not have been set ----------------------
[usr_par] = usr_par_init_default_parameters(usr_par);


%- get source and material ------------------------------------------------
noise_source = make_noise_source('no');
structure = define_material_parameters('no');


%- loop over reference stations -------------------------------------------
if (~exist(filename('correlations', n_ref), 'file'))
    correlations = zeros(n_ref, n_rec, nt);
else
    correlations = parload(filename('correlations', n_ref));
end

misfit = 0;
gradient = zeros(nx, nz, 2);

tic
for i_ref = 1:n_ref
    
    
    % each reference station will act as a source once --------------------
    src = usr_par.network.ref_stat(i_ref,:);
    rec = usr_par.network.array(~ismember(usr_par.network.array, src, 'rows'),:);
    
    
    % calculate Green function --------------------------------------------
    if (strcmp(usr_par.type, 'source'))
        
        if (~exist(filename('G_fft', i_ref), 'file'))
            if(usr_par.verbose); fprintf('ref %i: calculate Green function\n', i_ref); end
            G_fft = run1_forward_green(structure, src, 0);
            parsave(filename('G_fft', i_ref), G_fft)
        else
            if(usr_par.verbose); fprintf('ref %i: load pre-computed Green function\n', i_ref); end
            G_fft = parload(filename('G_fft', i_ref));
        end
        
    else
        if(usr_par.verbose); fprintf('ref %i: calculate Green function\n', i_ref); end
        [G_fft, G] = run1_forward_green(structure, src, 1);
    end
    
    
    % calculate correlation -----------------------------------------------
    if (strcmp(usr_par.type, 'source'))
        
        if (~exist(filename('correlations', n_ref), 'file'))
            if(usr_par.verbose); fprintf('ref %i: calculate correlations\n', i_ref); end
            correlations(i_ref,:,:) = run2_forward_correlation(structure, noise_source, G_fft, src, rec, 0);
        else
            if(usr_par.verbose); fprintf('ref %i: use pre-computed correlations\n', i_ref); end
        end
        
    else
        if(usr_par.verbose); fprintf('ref %i: calculate correlations\n', i_ref); end
        [correlations(i_ref,:,:), C] = run2_forward_correlation(structure, noise_source, G_fft, src, rec, 1);
    end
    
    
    % calculate misfits and adjoint source functions ----------------------
    [misfit_iref, adjstf_iref] = make_adjoint_sources( ...
        reshape(correlations(i_ref,:,:), [], nt), reshape(usr_par.data.correlations(i_ref,:,:), [], nt), ...
        t, usr_par.measurement.type, src, rec, usr_par.measurement.mode);
    
    
    % sum up misfits for all reference stations ---------------------------
    misfit = misfit + sum(misfit_iref);
    
    
    % calculate gradient and sum them up for all reference stations -------
    if (strcmp(usr_par.type, 'source'))
        
        if(usr_par.verbose); fprintf('ref %i: calculate source kernel\n', i_ref); end
        gradient_iref = run3_adjoint(structure, noise_source, G_fft, src, rec, adjstf_iref, [], 0);
        gradient = gradient + gradient_iref;
        
    else
        
        if(usr_par.verbose); fprintf('ref %i: calculate structure kernel - part 1\n', i_ref); end
        [gradient_iref_1, adjoint_state] = run3_adjoint(structure, noise_source, [], src, rec, adjstf_iref, C, 1);
        
        if(usr_par.verbose); fprintf('ref %i: calculate structure kernel - part 2\n', i_ref); end
        gradient_iref_2 = run3_adjoint(structure, noise_source, [], src, rec, adjoint_state, G, 0);
        
        gradient = gradient + gradient_iref_1 + gradient_iref_2;
        
    end
    
    
    if(usr_par.verbose); fprintf('ref %i: done\n', i_ref); end
    
    
end
toc


%- save correlations ------------------------------------------------------
if(usr_par.verbose); fprintf('misfit:   %15.10f\n', misfit); end
if (~exist(filename('correlations', n_ref), 'file'))
    save(filename('correlations', n_ref), 'correlations', 't')
end


%- plot kernel ------------------------------------------------------------
plot_kernel(gradient, usr_par)


