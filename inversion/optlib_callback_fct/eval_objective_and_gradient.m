
function [j, g, c_all] = eval_objective_and_gradient(m, ModRandString, usr_par, only_objective)
 
if( nargin < 4 )
    only_objective = false;
end

[~,~,nx,nz,~,~,~,model_type,~,n_basis_fct] = input_parameters();

%- redirect optimization variable x and initialize kernel structures
if( strcmp( usr_par.type, 'source') )
    
    % loading of spectra important for n_basis_fct=0, i.e. one map for each noise source
    source_dist = map_m_to_parameters(m, usr_par);
    [~, spectrum] = make_noise_source();
    
    % load structure that is assumed for the source inversion
    [mu, rho] = define_material_parameters(nx,nz,model_type);
    
elseif( strcmp( usr_par.type, 'structure') )
    
    % get source that is assumed for the structure inversion
    [source_dist, spectrum] = make_noise_source();
    
    % get mu from v, which is our optimization variable (relative parameterization)
    mu = map_m_to_parameters(m, usr_par);
    [~,rho] = define_material_parameters(nx,nz,model_type);
    
    % rho = map_m_to_parameters(m, usr_par);
    % [mu,~] = define_material_parameters(nx,nz,model_type);
    
end



%- loop over reference stations
j = 0;
t = usr_par.data.t;
nt = length(t);
n_ref = size( usr_par.network.ref_stat, 1 );
n_rec = size( usr_par.network.array, 1 ) - 1;

c_it = zeros( n_ref, n_rec, nt );

if(n_basis_fct == 0)
    n_basis_fct = 1;
end

grad_parameters = zeros( nx, nz, n_basis_fct);


parfor i = 1:n_ref
    
    
    % each reference station will act as a source once
    src = usr_par.network.ref_stat(i,:);
    rec = usr_par.network.array( find( ~ismember(usr_par.network.array, src, 'rows') ) ,:);
    
    
    % load or calculate Green function
    if( strcmp( usr_par.type, 'source') && exist(['../output/interferometry/G_2_' num2str(i) '.mat'], 'file') )
        
        G_2 = load_G_2( ['../output/interferometry/G_2_' num2str(i) '.mat'] );
        
    else
        
        if( strcmp( usr_par.type, 'source') )
            
            if( strcmp( usr_par.use_mex, 'yes') )
                [G_2] = run_forward1_green_mex(mu, rho, src, 0);
            else
                [G_2] = run_forward1_green(mu, rho, src, 0);
            end            
            parsave( ['../output/interferometry/G_2_' num2str(i) '.mat'], G_2 )

        else
            
            if( strcmp( usr_par.use_mex, 'yes') )
                [G_2, G_2_dxu_time, G_2_dzu_time] = run_forward1_green_mex(mu, rho, src, 1);
            else
                [G_2, G_2_dxu_time, G_2_dzu_time] = run_forward1_green(mu, rho, src, 1);
            end
        
            % parsave( ['../output/interferometry/G_2_' num2str(i) '.mat'], G_2 )
            % parsave( ['../output/interferometry/G_2_dxu_time_' num2str(i) '.mat'], G_2_dxu_time )
            % parsave( ['../output/interferometry/G_2_dzu_time_' num2str(i) '.mat'], G_2_dzu_time )
            
        end
    end
    
    
    % calculate correlation
    if( strcmp( usr_par.type, 'source') )
        
        if( strcmp( usr_par.use_mex, 'yes') )
            [c_it(i,:,:)] = run_forward2_correlation_mex( mu, rho, G_2, spectrum, source_dist, rec, 0, usr_par.debug.df );
        else
            [c_it(i,:,:)] = run_forward2_correlation( mu, rho, G_2, spectrum, source_dist, rec, 0, usr_par.debug.df );
        end
        
    elseif( strcmp( usr_par.type, 'structure') )
        
        if( strcmp( usr_par.use_mex, 'yes') )
            [c_it(i,:,:), ~, C_2_dxu_time, C_2_dzu_time] = run_forward2_correlation_mex( mu, rho, G_2, spectrum, source_dist, rec, 1, usr_par.debug.df );
        else
            [c_it(i,:,:), ~, C_2_dxu_time, C_2_dzu_time] = run_forward2_correlation( mu, rho, G_2, spectrum, source_dist, rec, 1, usr_par.debug.df );
        end
        
    end
    
    
    % filter correlations if wanted
    if( strcmp( usr_par.filter.apply_filter, 'yes') )
        c_data_iref = filter_correlations( usr_par.data.c_data( (i-1)*n_rec + 1 : i*n_rec, : ), t, usr_par.filter.f_min, usr_par.filter.f_max, 4 );
        c_it(i,:,:) = filter_correlations( reshape(c_it(i,:,:),[],nt), t, usr_par.filter.f_min, usr_par.filter.f_max, 4 );
               
    else
        c_data_iref = usr_par.data.c_data( (i-1)*n_rec + 1 : i*n_rec, : );
        
    end
    
    
    % calculate misfit and adjoint source function
    [j_n, adstf] = make_adjoint_sources_inversion( reshape(c_it(i,:,:),[],nt), c_data_iref, t, usr_par.veldis, usr_par.measurement, src, rec );

    
    if( only_objective == false )
        
        % calculate source kernel
        if( strcmp( usr_par.type, 'source') )
            
            if( strcmp( usr_par.use_mex, 'yes') )
                grad_parameters_i = run_noise_source_kernel_mex( mu, rho, G_2, spectrum, adstf, rec );
            else
                grad_parameters_i = run_noise_source_kernel( mu, rho, G_2, spectrum, adstf, rec );
            end
            
        % calculate structure kernel
        elseif( strcmp( usr_par.type, 'structure') )
            
            if( strcmp( usr_par.use_mex, 'yes') )
                
                % first run
                [grad_mu_i_1, adjoint_state_1] = run_noise_structure_kernel_mex( mu, rho, C_2_dxu_time, C_2_dzu_time, complex(adstf), rec, spectrum, source_dist );
                
                % second run
                [grad_mu_i_2] = run_noise_structure_kernel_mex( mu, rho, G_2_dxu_time, G_2_dzu_time, adjoint_state_1, rec, spectrum, source_dist );
                
            else
                
                % first run
                [grad_mu_i_1, adjoint_state_1] = run_noise_structure_kernel( mu, rho, C_2_dxu_time, C_2_dzu_time, complex(adstf), rec, spectrum, source_dist );
                
                % second run
                [grad_mu_i_2] = run_noise_structure_kernel( mu, rho, G_2_dxu_time, G_2_dzu_time, adjoint_state_1, rec, spectrum, source_dist );
                
            end
            
            % sum up both contributions
            grad_parameters_i = grad_mu_i_1 + grad_mu_i_2;
            
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


%- kernel treatment, i.e. clipping/smoothing
grad_parameters = treat_kernel( grad_parameters, usr_par );


%- map grad_parameters to grad_m
g = map_gradparameters_to_gradm( m, grad_parameters, usr_par );


end


