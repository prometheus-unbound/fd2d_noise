
function [misfit, grad_m, c_iteration] = eval_objective_and_gradient( m, ModRandString, usr_par )


[~, ~, nx, nz] = input_parameters();


% inversion toolbox requires m to be a vector
m_parameters = reshape( m, nx, nz, [] );


% material parameters
structure.mu = m_parameters(:,:,1);
structure.rho = usr_par.structure.rho;


% get source and material
noise_source.distribution = m_parameters(:,:,2);
noise_source.spectrum = usr_par.noise_source.spectrum;


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
        c_iref(i_ref,:,:) = run2_forward_correlation_mex( structure, noise_source, G_fft, src, rec, 0 );
    else
        [c_iref(i_ref,:,:), C] = run2_forward_correlation_mex( structure, noise_source, G_fft, src, rec, 1 );
    end
    
    
    % calculate misfits and adjoint source functions
    [misfit_iref, adjstf_iref] = make_adjoint_sources( reshape(c_iref(i_ref,:,:),[],nt), usr_par.data.correlations( (i_ref-1)*n_rec + 1 : i_ref*n_rec, : ), t, usr_par.measurement.type, src, rec, usr_par.measurement.mode );
    
    
    % sum up misfits for all reference stations
    misfit = misfit + sum(misfit_iref);
    
    
    % calculate kernel and sum them up for all reference stations
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


% set gradient in absorbing boundaries to zero
absbound = init_absbound();
for i = 1:2
   gradient(:,:,i) = double( absbound == 1 ) .* gradient(:,:,i); 
end


% inversion toolbox requires a vector
grad_m = reshape( gradient, [], 1 );


end

