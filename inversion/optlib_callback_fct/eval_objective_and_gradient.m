
function [j, g, c_all] = eval_objective_and_gradient(m, ModRandString, usr_par, only_objective)
 

if( nargin < 4 )
    only_objective = false;
end


[~,~,nx,nz,~,~,~,model_type,~,n_basis_fct] = input_parameters();
if(n_basis_fct == 0)
    n_basis_fct = 1;
end


%- redirect optimization variable x and initialize kernel structures
parameters = map_m_to_parameters(m, usr_par);

% loading of spectra important for n_basis_fct=0, i.e. one map for each noise source
source_dist = parameters(:,:,1:n_basis_fct);
[~, spectrum] = make_noise_source();

% material parameters
mu = parameters(:,:,end);
[~,rho] = define_material_parameters(nx,nz,model_type);


%- loop over reference stations
j = 0;
t = usr_par.data.t;
nt = length(t);
n_ref = size( usr_par.network.ref_stat, 1 );
n_rec = size( usr_par.network.array, 1 ) - 1;
c_it = zeros( n_ref, n_rec, nt );
grad_parameters = zeros( nx, nz, n_basis_fct + 1 );


parfor i = 1:n_ref
    
    
    % each reference station will act as a source once
    src = usr_par.network.ref_stat(i,:);
    rec = usr_par.network.array( find( ~ismember(usr_par.network.array, src, 'rows') ) ,:);
    
    
    % calculate Green function 
    if( strcmp( usr_par.use_mex, 'yes') )
        [G_2, G_2_dxu_time, G_2_dzu_time] = run_forward1_green_mex(mu, rho, src, 1);
    else
        [G_2, G_2_dxu_time, G_2_dzu_time] = run_forward1_green(mu, rho, src, 1);
    end
            
    
    % calculate correlation        
    if( strcmp( usr_par.use_mex, 'yes') )
        [c_it(i,:,:), ~, C_2_dxu_time, C_2_dzu_time] = run_forward2_correlation_mex( mu, rho, G_2, spectrum, source_dist, rec, 1, usr_par.debug.df );
    else
        [c_it(i,:,:), ~, C_2_dxu_time, C_2_dzu_time] = run_forward2_correlation( mu, rho, G_2, spectrum, source_dist, rec, 1, usr_par.debug.df );
    end
    
    
    % filter correlations if wanted
    if( strcmp( usr_par.filter.apply_filter, 'yes') )
        c_data_iref = filter_correlations( usr_par.data.c_data( (i-1)*n_rec + 1 : i*n_rec, : ), t, usr_par.filter.f_min, usr_par.filter.f_max, 4 );
        c_it(i,:,:) = filter_correlations( reshape(c_it(i,:,:),[],nt), t, usr_par.filter.f_min, usr_par.filter.f_max, 4 );
               
    else
        c_data_iref = usr_par.data.c_data( (i-1)*n_rec + 1 : i*n_rec, : );
        
    end
    
    
    % calculate misfits and adjoint source functions
    [j_n_source, adstf_source] = make_adjoint_sources_inversion( reshape(c_it(i,:,:),[],nt), c_data_iref, t, usr_par.veldis, usr_par.measurement.source, src, rec );
    [j_n_structure, adstf_structure] = make_adjoint_sources_inversion( reshape(c_it(i,:,:),[],nt), c_data_iref, t, usr_par.veldis, usr_par.measurement.structure, src, rec );
    
    
    % build up total misfit and adjoint source time function
    j_n = (1.0 - usr_par.kernel.weighting) * j_n_source + usr_par.kernel.weighting * j_n_structure;
    adstf = (1.0 - usr_par.kernel.weighting) * adstf_source +  usr_par.kernel.weighting * adstf_structure;
    
    
    if( only_objective == false )
        
        % calculate kernel kernel
        if( strcmp( usr_par.use_mex, 'yes') )
            
            % first run
            [grad_i_1, adjoint_state_1] = run_noise_adjoint_mex( mu, rho, C_2_dxu_time, C_2_dzu_time, complex(adstf), rec, spectrum, source_dist, G_2 );
            
            % second run
            [grad_i_2] = run_noise_adjoint_mex( mu, rho, G_2_dxu_time, G_2_dzu_time, adjoint_state_1, rec, spectrum, source_dist, G_2 );
            
        else
            
            % first run
            [grad_i_1, adjoint_state_1] = run_noise_adjoint( mu, rho, C_2_dxu_time, C_2_dzu_time, complex(adstf), rec, spectrum, source_dist, G_2 );
            
            % second run
            [grad_i_2] = run_noise_adjoint( mu, rho, G_2_dxu_time, G_2_dzu_time, adjoint_state_1, rec, spectrum, source_dist, G_2 );
            
        end
        
        % sum up both contributions
        grad_parameters_i = grad_i_1 + grad_i_2;
        
        
        % sum up kernels
        grad_parameters = grad_parameters + grad_parameters_i;
        
    end
    
    
    % sum up misfits
    j = j + sum(j_n);
    
    
end


fprintf('misfit: %f\n',j)

%%% TEST REGULARIZATON
j = j + usr_par.regularization.alpha/2 * sum( usr_par.regularization.weighting .* (m - usr_par.m0).^2 ); 
fprintf('misfit with regularization: %f\n',j)


%- reorganize correlation vector
c_all = zeros( n_ref*n_rec, length(t) );
for i = 1:n_ref
    c_all( (i-1)*n_rec + 1 : i*n_rec, :) = c_it(i,:,:);
end


if( only_objective == true )
    g = 0;
    return    
end


%- kernel treatment, i.e. clipping
grad_parameters = treat_kernel( grad_parameters, usr_par );


%- map grad_parameters to grad_m
g = map_gradparameters_to_gradm( m, grad_parameters, usr_par );


%%% TEST REGULARIZATON
g = g + usr_par.regularization.alpha * usr_par.regularization.weighting .* ( m - usr_par.m0 );


end


