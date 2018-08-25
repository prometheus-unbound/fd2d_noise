
function [j, g, c_all] = eval_objective_and_gradient( m, ModRandString, usr_par, only_objective )
 

if( nargin < 4 )
    only_objective = false;
end


%- get configuration
[~,~,~,~,~,~,~,model_type, source_type, n_basis_fct] = input_parameters();

%- redirect optimization variable x and initialize kernel structures
m_parameters = map_m_to_parameters(m, usr_par);

% loading of spectra important for n_basis_fct=0, i.e. one map for each noise source
source_dist = m_parameters(:,:,1:end-1);
[~, spectrum] = make_noise_source(source_type, n_basis_fct);

% material parameters
mu = m_parameters(:,:,end);
[~,rho] = define_material_parameters( usr_par.config.nx, usr_par.config.nz, model_type );


%- loop over reference stations
j = 0; j_source = 0; j_structure = 0;
t = usr_par.data.t;
nt = length(t);
n_ref = size( usr_par.network.ref_stat, 1 );
n_rec = size( usr_par.network.array, 1 ) - 1;
c_it = zeros( n_ref, n_rec, nt );
grad_parameters = 0.0 * m_parameters;


parfor i = 1:n_ref
    
    
    % each reference station will act as a source once
    src = usr_par.network.ref_stat(i,:);
    rec = usr_par.network.array( find( ~ismember(usr_par.network.array, src, 'rows') ) ,:);
    
    
    % calculate Green function 
    if( strcmp( usr_par.type, 'source') && exist(['../output/interferometry/G_fft_' num2str(i) '.mat'], 'file') )
        G_fft = parload( ['../output/interferometry/G_fft_' num2str(i) '.mat'] );
    else
        
        if( strcmp( usr_par.type, 'source' ) )
            [G_fft] = run_forward1_green_mex(mu, rho, src, [], 0, [], single([]) );
            parsave( ['../output/interferometry/G_fft_' num2str(i) '.mat'], G_fft )
        else
            [G_fft, G] = run_forward1_green_mex(mu, rho, src, [], 1, [], single([]));
        end
        
    end
            
    
    % calculate correlation
    if( strcmp( usr_par.type, 'source' ) )
        [c_it(i,:,:)] = run_forward2_correlation_mex( mu, rho, G_fft, spectrum, source_dist, rec, 0, [], single([]) );
    else
        [c_it(i,:,:), C] = run_forward2_correlation_mex( mu, rho, G_fft, spectrum, source_dist, rec, 1, [], single([]) );
    end
    
    
    % filter correlations if wanted
    if( strcmp( usr_par.filter.apply_filter, 'yes' ) )
        c_data_iref = filter_correlations( usr_par.data.c_data( (i-1)*n_rec + 1 : i*n_rec, : ), t, usr_par.filter.f_min, usr_par.filter.f_max, 4 );
        c_it(i,:,:) = filter_correlations( reshape(c_it(i,:,:),[],nt), t, usr_par.filter.f_min, usr_par.filter.f_max, 4 );       
    else
        c_data_iref = usr_par.data.c_data( (i-1)*n_rec + 1 : i*n_rec, : );
    end
    
    
    % calculate misfits and adjoint source functions
    [j_n_source, adstf_source] = make_adjoint_sources( reshape(c_it(i,:,:),[],nt), c_data_iref, 0*c_data_iref, t, usr_par.veldis, usr_par.measurement.source, src, rec, '1st' );
    [j_n_structure, adstf_structure] = make_adjoint_sources( reshape(c_it(i,:,:),[],nt), c_data_iref, 0*c_data_iref, t, usr_par.veldis, usr_par.measurement.structure, src, rec, '1st' );
    
    
    % build up total misfit and adjoint source time function
    j_n = (1.0 - usr_par.kernel.weighting) * j_n_source + usr_par.kernel.weighting * j_n_structure;
    adstf = (1.0 - usr_par.kernel.weighting) * adstf_source +  usr_par.kernel.weighting * adstf_structure;
    
    
    if( only_objective == false )
        
        % calculate kernel kernel
        if( strcmp( usr_par.type, 'source' ) )
            
            % mode == 0, do not need Fourier transformed adjoint state
            [grad_i_1] = run_noise_adjoint_mex( mu, rho, single([]), complex(adstf), rec, [], spectrum, [], G_fft, 0, [], single([]) );
            [grad_i_2] = 0.0 * grad_i_1;
            
        elseif( strcmp( usr_par.type, 'structure' ) )
            
            % mode == 1, need Fourier transformed adjoint state
            [grad_i_1, adjoint_state_1] = run_noise_adjoint_mex( mu, rho, C, complex(adstf), rec, [], [], [], [], 1, [], single([]) );
            [grad_i_2] = run_noise_adjoint_mex( mu, rho, G, adjoint_state_1, rec, [], spectrum, source_dist, [], 0, [], single([]) );
            
        else
            
            % mode == 1, need Fourier transformed adjoint state
            [grad_i_1, adjoint_state_1] = run_noise_adjoint_mex( mu, rho, C, complex(adstf), rec, [], spectrum, [], G_fft, 1, [], single([]) );
            [grad_i_2] = run_noise_adjoint_mex( mu, rho, G, adjoint_state_1, rec, [], spectrum, source_dist, [], 0, [], single([]) );
            
        end      
        
        
        % sum up kernels
        grad_parameters = grad_parameters + grad_i_1 + grad_i_2;
        
        
    end
    
    
    % sum up misfits
    j = j + sum(j_n);
    j_source = j_source + sum(j_n_source);
    j_structure = j_structure + sum(j_n_structure);
    
    
end


fprintf( 'misfit total:      %15.10f\n', j )
fprintf( 'misfit source:     %15.10f\n', j_source )
fprintf( 'misfit structure:  %15.10f\n', j_structure )


%- reorganize correlation vector
c_all = zeros( n_ref*n_rec, length(t) );
for i = 1:n_ref
    c_all( (i-1)*n_rec + 1 : i*n_rec, :) = c_it(i,:,:);
end


if( only_objective == true )
    
    j = regularization( j, 0, m, usr_par );
    fprintf('misfit with regu.: %15.10f\n',j)
    g = 0;    
    return    
    
end


%- map grad_parameters to grad_m
norm_first = sum(sum(sum(grad_parameters)))
g = map_gradparameters_to_gradm( m, grad_parameters, usr_par );

%- add regularization
[ j, g ] = regularization( j, g, m, usr_par );
fprintf('misfit with regu.: %15.10f\n',j)


end

