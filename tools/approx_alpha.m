
function [ alpha, beta ] = approx_alpha( misfit )

    [~,~,nx,nz,~,~,~,~,~,n_basis_fct] = input_parameters();
    [usr_par] = usr_par_init_default_parameters_lbfgs([]);
    
    
    if(n_basis_fct == 0)
        m_parameters = zeros( nx, nz, 2 );
        m_parameters_0 = ones( nx, nz, 2 );
    else
        m_parameters = zeros( nx, nz, n_basis_fct+1 );
        m_parameters_0 = zeros( nx, nz, n_basis_fct+1 );
    end
    
    
    m_parameters(:,:,1:end-1) = make_noise_source( 'gaussian', n_basis_fct );
    m_parameters(:,:,end) = define_material_parameters( nx, nz, 200 );
    
    m_parameters_0(:,:,1:end-1) = make_noise_source( 'homogeneous', n_basis_fct );
    m_parameters_0(:,:,end) = define_material_parameters( nx, nz, 1 );


    % convert to optimization variable
    m = map_parameters_to_m( m_parameters, usr_par );
    m_0 = map_parameters_to_m( m_parameters_0, usr_par );
    
    
    % approximate alpha and beta
    if(n_basis_fct == 0)
        alpha = 2 * 0.01 * misfit / sum( weighting( nx, nz ) .* ( m( 1 : end - nx*nz, 1) - m_0( 1 : end - nx*nz, 1) ).^2 );
    else
        alpha = 2 * 0.01 * misfit / sum( repmat(weighting( nx, nz ), n_basis_fct, 1) .* ( m( 1 : end - nx*nz, 1) - m_0( 1 : end - nx*nz, 1) ).^2 );
    end
    
    beta  = 2 * 0.01 * misfit / sum( weighting( nx, nz ) .* ( m(end - nx*nz + 1 : end , 1) - m_0(end - nx*nz + 1  : end , 1 ) ).^2 );

end

