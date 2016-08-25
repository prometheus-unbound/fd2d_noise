
function [misfit, grad_m] = eval_objective_and_gradient(m, ModRandString, usr_par)


%- inversion toolbox requires m to be a vector ----------------------------
m_parameters = map_m_to_parameters( m, usr_par );


%- material parameters ----------------------------------------------------
structure.mu = m_parameters(:,:,1);
structure.rho = usr_par.structure.rho;


%- get source and material ------------------------------------------------
noise_source.distribution = m_parameters(:,:,2);
noise_source.spectrum = usr_par.noise_source.spectrum;


%- loop over reference stations -------------------------------------------
n_ref = size(usr_par.network.ref_stat, 1);
n_rec = size(usr_par.network.array, 1) - 1;
t = usr_par.data.t;
nt = length(t);
correlations = zeros(n_ref, n_rec, nt);
misfit = 0;
gradient = zeros(size(m_parameters, 1), size(m_parameters, 2), 2);


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
        if(usr_par.verbose); fprintf('ref %i: calculate correlations\n', i_ref); end
        correlations(i_ref,:,:) = run2_forward_correlation(structure, noise_source, G_fft, src, rec, 0);
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
    
    
    % completed calculations for this reference station -------------------
    if(usr_par.verbose); fprintf('ref %i: done\n', i_ref); end
    
    
end


%- print total misfit -----------------------------------------------------
if(usr_par.verbose); fprintf('misfit:   %15.10f\n', misfit); end


%- set gradient in absorbing boundaries to zero ---------------------------
absbound = init_absbound();
for i = 1:2
    gradient(:,:,i) = double(absbound == 1) .* gradient(:,:,i);
end


%- convert gradient w.r.t. parameters to gradient w.r.t. opt. variables ---
grad_m = map_gradparameters_to_gradm( m, gradient, usr_par );


end


