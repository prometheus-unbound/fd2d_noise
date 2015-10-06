
function [j, g, c_all] = eval_objective_and_gradient(m, ModRandString, usr_par, only_objective)
 
if( nargin < 4 )
    only_objective = false;
end

[~,~,nx,nz,~,~,~,model_type] = input_parameters();

%- redirect optimization variable x and initialize kernel structures
if( strcmp( usr_par.type, 'source') )
    
    % loading of spectra important for n_basis_fct=0, i.e. one map for each noise source
    source_dist = map_m_to_parameters(m, usr_par);
    [~,spectrum] = make_noise_source();
    
    % load structure that is assumed for the source inversion
    mu = define_material_parameters(nx,nz,model_type);
    
elseif( strcmp( usr_par.type, 'structure') )
    
    % get source that is assumed for the structure inversion
    [source_dist,spectrum] = make_noise_source();
    
    % get mu from v, which is our optimization variable (relative parameterization)
    mu = map_m_to_parameters(m, usr_par);
    
end



%- loop over reference stations
j = 0;
t = usr_par.data.t;
nt = length(t);
n_ref = size( usr_par.network.ref_stat, 1 );
n_rec = size( usr_par.network.array, 1 ) - 1;

c_it = zeros( n_ref, n_rec, nt );
grad_parameters = zeros( nx, nz );

for i = 1:n_ref
    
    
    % each reference station will act as a source once
    src = usr_par.network.ref_stat(i,:);
    rec = usr_par.network.array( find( ~ismember(usr_par.network.array, src, 'rows') ) ,:);
    
    
    % load or calculate Green function
    if( strcmp( usr_par.type, 'source') && exist(['../output/interferometry/G_2_' num2str(i) '.mat'], 'file') )
        G_2 = load_G_2( ['../output/interferometry/G_2_' num2str(i) '.mat'] );
        
    else
        [G_2] = run_forward_green_fast_mex(mu, src);
        
        if( strcmp( usr_par.type, 'source') )
            parsave( ['../output/interferometry/G_2_' num2str(i) '.mat'], G_2 )
        end
    end
    
    
    % calculate correlation
    if( strcmp( usr_par.type, 'source') )
        [c_it(i,:,:)] = run_forward_correlation_fast_mex( G_2, source_dist, spectrum, mu, rec, 0, usr_par.debug.df );
        
    elseif( strcmp( usr_par.type, 'structure') )
        [c_it(i,:,:), ~, C_2_dxv, C_2_dzv] = run_forward_correlation_fast_mex( G_2, source_dist, spectrum, mu, rec, 1, usr_par.debug.df );
        
    end
    
    
    % filter correlations if wanted
    if( strcmp( usr_par.filter.apply_filter, 'yes') )
        c_data_iref = filter_correlations( usr_par.data.c_data( (i-1)*n_rec + 1 : i*n_rec, : ), t, usr_par.filter.f_min, usr_par.filter.f_max );
        c_it(i,:,:) = filter_correlations( reshape(c_it(i,:,:),[],nt), t, usr_par.filter.f_min, usr_par.filter.f_max );
        
    else
        c_data_iref = usr_par.data.c_data( (i-1)*n_rec + 1 : i*n_rec, : );
        
    end
    
    
    % calculate misfit and adjoint source function
    [j_n, adstf] = make_adjoint_sources_inversion( reshape(c_it(i,:,:),[],nt), c_data_iref, t, usr_par.veldis, usr_par.measurement, src, rec );

    
    if( only_objective == false )
        
        % calculate source kernel
        if( strcmp( usr_par.type, 'source') )

            grad_parameters_i = run_noise_source_kernel_fast_mex( G_2, mu, spectrum, adstf, rec );
            
        % calculate structure kernel
        elseif( strcmp( usr_par.type, 'structure') )
            
            grad_parameters_i = run_noise_mu_kernel_fast_mex( C_2_dxv, C_2_dzv, mu, adstf, rec );
            
        end                
        
        % sum up kernels
        grad_parameters = grad_parameters + grad_parameters_i;
        
    end
    
    
    % sum up misfits
    j = j + sum(j_n);
    
    
end


fprintf('misfit: %f\n',j)


%- reorganize correlation vector
c_all = zeros( n_ref*n_rec, length(t) );
for i = 1:n_ref
    c_all( (i-1)*n_rec + 1 : i*n_rec, :) = c_it(i,:,:);
end


if( only_objective == true )
    g = 0;
    return    
end


%- sum source kernel (all frequencies or in frequency bands)
if( strcmp( usr_par.type, 'source') )
    grad_parameters = sum_source_kernel( grad_parameters );
end


%- kernel treatment, i.e. clipping/smoothing
grad_parameters = treat_kernel( grad_parameters, usr_par );


%- map grad_parameters to grad_m
g = map_gradparameters_to_gradm( m, grad_parameters, usr_par );


end


