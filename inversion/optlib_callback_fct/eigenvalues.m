
function [Hdm] = eigenvalues( dm )
    

    global counter
    counter = counter + 1;
    fprintf('%i\n',counter)

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % user input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    usr_par.network = load( '../output/interferometry/array_16_ref_small.mat' );
    usr_par.data = load( '../output/interferometry/data_16_ref_0_small_gaussian.mat' );
    
    % usr_par.kernel.imfilter.source = fspecial('gaussian',[40 40], 20);
    usr_par.kernel.sigma.source = [1 1];
    
    usr_par.use_mex = 'yes';
    usr_par.type = 'source';
    usr_par.kernel.weighting = 0.0;
    usr_par.measurement.source = 'waveform_difference';
    usr_par.measurement.structure = 'waveform_difference';
    usr_par.verbose = 'yes';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate Hessian vector product
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [usr_par] = usr_par_init_default_parameters_lbfgs( usr_par );
    [ ~, ~, nx, nz, ~, ~, ~, model_type, source_type, n_basis_fct ] = input_parameters();
    

    % set up initial model
    if( n_basis_fct == 0)
        m_parameters = zeros( nx, nz, 2 );
        m_parameters(:,:,1) = make_noise_source( source_type, n_basis_fct );
    else
        m_parameters = zeros( nx, nz, n_basis_fct + 1 );
        m_parameters(:,:,1:n_basis_fct) = make_noise_source( source_type, n_basis_fct );
    end
    m_parameters(:,:,end) = define_material_parameters( nx, nz, model_type );


    % convert to optimization variable
    m = map_parameters_to_m( m_parameters, usr_par );
    
    
    % set test vectors in absorbing boundary region to zero
    [absbound] = reshape( init_absbound() , [], 1 );
    for i = 1:size(dm)/(nx*nz)
       dm( ((i-1)*nx*nz+1):i*nx*nz, 1 ) = double( absbound == 1 ) .* dm( ((i-1)*nx*nz+1):i*nx*nz, 1 ); 
    end
    
    
    % set dm to zero for source and structure case
    if( numel(m_parameters) == numel(dm) )
        
        if( strcmp( usr_par.type, 'source' ) )
            dm( end-nx*nz+1:end, 1 ) = 0 * dm( end-nx*nz+1:end, 1 );
        elseif( strcmp( usr_par.type, 'structure' ) )
            dm( 1:end-nx*nz, 1 ) = 0 * dm( 1:end-nx*nz, 1 );
        end
        
    end

    
    % calculate Hessian vector product
    if( numel(m_parameters) == numel(dm) && ~isempty(find( dm( end-nx*nz+1:end, 1 ), 1 ))  &&  ~isempty(find( dm( 1:end-nx*nz, 1 ), 1 ))  )
        
        dm_tmp = 0 * dm;
        dm_tmp( 1:end-nx*nz, 1 ) = dm( 1:end-nx*nz, 1 );
        Hdm = eval_hessian_vector_product( m, dm_tmp, optlib_generate_random_string(8), usr_par );
        
        dm_tmp = 0 * dm;
        dm_tmp( end-nx*nz+1:end, 1 ) = dm( end-nx*nz+1:end, 1 );
        Hdm_tmp = eval_hessian_vector_product( m, dm_tmp, optlib_generate_random_string(8), usr_par );
        
        Hdm = Hdm + Hdm_tmp;
        
    else
        
        % Hdm = eval_hessian_vector_product( m, dm, optlib_generate_random_string(8), usr_par );
        
        tic
        Hdm = eval_hessian_vector_product_source( m, dm, optlib_generate_random_string(8), usr_par );
        toc
        
    end
    
    
    % only return part that is of interest
    if( strcmp( usr_par.type, 'source' ) )
        Hdm = Hdm( 1:end-nx*nz, 1 );
    elseif( strcmp( usr_par.type, 'structure') )
        Hdm = Hdm( end-nx*nz+1:end, 1 );        
    end


end
