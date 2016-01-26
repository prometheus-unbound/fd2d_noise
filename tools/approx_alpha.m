
function [ alpha ] = approx_alpha( misfit )

    [~,~,nx,nz,~,~,~,model_type,~,n_basis_fct] = input_parameters();
    [usr_par] = usr_par_init_default_parameters_lbfgs([]);
    
    % set up initial model
    if( n_basis_fct == 0)
        m_parameters = zeros(nx, nz, 2);
        m_parameters(:,:,1) = make_noise_source();
        
        m_parameters_0 = zeros(nx, nz, 2);
        m_parameters_0(:,:,1) = ones(nx,nz);
    else
        m_parameters = zeros(nx, nz, n_basis_fct+1);
        m_parameters(:,:,1:n_basis_fct) = make_noise_source();
        
        m_parameters_0 = zeros(nx, nz, n_basis_fct+1);
        m_parameters_0(:,:,1:n_basis_fct) = ones(nx,nz,n_basis_fct);
    end
    
    m_parameters(:,:,end) = define_material_parameters(nx,nz,9999);
    m_parameters_0(:,:,end) = define_material_parameters(nx,nz,1);


    % convert to optimization variable
    m = map_parameters_to_m(m_parameters,usr_par);
    m_0 = map_parameters_to_m(m_parameters_0,usr_par);
    
    
    % approximate alpha
    alpha = 2 * 0.01 * misfit / sum( weighting( nx, nz, n_basis_fct ) .* (m - m_0).^2 );

end

