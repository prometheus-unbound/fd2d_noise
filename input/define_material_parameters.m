function [mu, rho] = define_material_parameters( nx, nz, model_type, make_plots )

%==========================================================================
% generate material parameters mu [N/m^2] and rho [kg/m^3]
%
% input: grid points of the velocity and density field in x-direction (nx) and z-direction (nz)
%        model_type
%        for "stand-alone" it is convenient to use make_plots as 'yes'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define make_plots if not specified
if( nargin < 4 )
    make_plots = 'no';
end


if (model_type==1)
    
    rho = 3000.0*ones(nx,nz);
    mu = 4.8e10*ones(nx,nz);
    
    
elseif (model_type==2)
    
    rho = 3000.0*ones(nx,nz);
    mu = 4.8e10*ones(nx,nz);
    
    % rho(98:102,123:127) = rho(98:102,123:127) + 2000.0;
    % mu(140:160,140:160) = mu(140:160,140:160) + 2.0e10;
    
    [Lx,Lz,nx,nz] = input_parameters();
    [X,Z] = define_computational_domain(Lx,Lz,nx,nz);
    x_sourcem = Lx/2.0;
    z_sourcem = Lz/2.0;
    x_width = Lx/10;
    
    mu = mu - 8.0e9 * ( exp( -( (X-x_sourcem).^2 + (Z-z_sourcem).^2 ) / (x_width)^2 ) );
    
    
elseif (model_type==3)
    
    rho = 3000.0*ones(nx,nz);
    mu = 4.8e10*ones(nx,nz);
    
    % rho(round(nx/2*1.25):end,:) = rho(round(nx/2*1.25):end,:) + 200.0;
    mu(round(nx/2*1.25):end,:) = mu(round(nx/2*1.25):end,:) + 0.8e10;
    
    
elseif (model_type==4)
    
    rho = 3000.0*ones(nx,nz);
    mu = 2.8e10*ones(nx,nz);
    
    rho(1:round(nx/2),:) = rho(1:round(nx/2),:) + 200.0;
    mu(1:round(nx/2),:) = mu(1:round(nx/2),:) + 3.5e10;
    
    rho(98:102,123:127)=rho(98:102,123:127)+2000.0;
    
    
elseif (model_type==5)
    
    rho = 3000.0*ones(nx,nz);
    mu = 2.0e10*ones(nx,nz);
    
    for k = 100:150
        mu(:,k) = mu(:,k) + (k-100)*10.0e8;
    end
    for k = 151:nz
        mu(:,k) = mu(:,150);
    end
   
    
elseif (model_type==6)
    
    rho = 3000.0*ones(nx,nz);
    mu = 2.0e10*ones(nx,nz);
    
    for k = 100:150
        mu(:,k) = mu(:,k) + (k-100)*10.0e8;
    end
    for k = 151:nz
        mu(:,k) = mu(:,150);
    end
    
    rho(98:102,123:127) = rho(98:102,123:127) + 2000.0;
    
    
elseif (model_type==7)
    
    rho = 3000.0*ones(nx,nz);
    mu = ones(nx,nz);
    mu(1:330,:) = 3.675e10;
    mu(331:end,:) = 2.7e10;

    
elseif ( model_type==666 || model_type==888 )
    
    rho = 3000.0 * ones(nx,nz);
    mu = 4.8e10 * ones(nx,nz);
    
    
elseif (model_type==999)
    
    rho = 3000.0*ones(nx,nz);
    mu = 4.8e10*ones(nx,nz);
    
    x_sourcem = [1.15e6 1.45e6];
    z_sourcem = [1.0e6 1.0e6];
    x_width = [1.1e5 1.1e5];
    z_width = [2.5e5 2.5e5];
    
    [Lx,Lz] = input_parameters();
    [X,Z] = define_computational_domain(Lx,Lz,nx,nz);
    
    % structure_1: 4.0e9
    % structure_2: 2.0e9
    % structure_3: 1.0e9
    for i=1:size(x_sourcem,2)
        mu = mu + (-1)^i * 1.0e10 * exp( -( (X-x_sourcem(i)).^2 ) / x_width(i)^2 )' .* exp( -( (Z-z_sourcem(i)).^2 ) / z_width(i)^2 )' ;
    end
    
    
elseif (model_type==9999)
    
    rho = 3000.0*ones(nx,nz);
    mu = 4.8e10*ones(nx,nz);
    
    x_sourcem = [2.3e5 2.9e5];
    z_sourcem = [2.0e5 2.0e5];
    x_width = [2.2e4 2.2e4];
    z_width = [5e4 5e4];
    
    [Lx,Lz] = input_parameters();
    [X,Z] = define_computational_domain(Lx,Lz,nx,nz);

    for i=1:size(x_sourcem,2)
        mu = mu + (-1)^i * 4.0e9 * exp( -( (X-x_sourcem(i)).^2 ) / x_width(i)^2 )' .* exp( -( (Z-z_sourcem(i)).^2 ) / z_width(i)^2 )' ;
    end

    
    
elseif (model_type==200)
    
    rho = 3000.0*ones(nx,nz);
    mu = 4.8e10*ones(nx,nz);
        
    x_sourcem = [1.275e6];
    z_sourcem = [0.975e6];
    x_width = [0.7e2];
    z_width = [0.7e2];
    
    % point1: 0.1e5, 0.1e5, 2e10
    % point2: 0.7e2, 0.7e2, 3.0e31
    
    [Lx,Lz] = input_parameters();
    [X,Z] = define_computational_domain(Lx,Lz,nx,nz);
    
    for i=1:size(x_sourcem,2)
        mu = mu + (-1)^i * 3.0e31 * exp( -( (X-x_sourcem(i)).^2 ) / x_width(i)^2 )' .* exp( -( (Z-z_sourcem(i)).^2 ) / z_width(i)^2 )' ;
    end
    
    
else

    load(['models/mu_' str(model_type)]);
    load(['models/rho_' str(model_type)]);
    
    
end



if( strcmp(make_plots,'yes') )
    
    [Lx,Lz,nx,nz,~,~,~,~,source_type,n_basis_fct] = input_parameters();   
    [X,Z] = define_computational_domain(Lx,Lz,nx,nz);
    
    load ../output/interferometry/array_16_ref.mat
    min_x = min(array(:,1));
    min_z = min(array(:,2));
    max_x = max(array(:,1));
    max_z = max(array(:,2));
    
    buffer = 1e5;
    pattern = double( X > (min_x-buffer) & X < (max_x+buffer) ) .* double( Z > (min_z-buffer) & Z < (max_z+buffer) );
    
    % array = [];
    
    if( n_basis_fct == 0 )
        m_parameters = zeros( nx, nz, 2);
    else
        m_parameters = zeros( nx, nz, n_basis_fct+1 );
    end
    
    m_parameters(:,:,1:end-1) = make_noise_source( source_type, n_basis_fct );
    m_parameters(:,:,end) = mu;
    
    
    usr_par.network = []; usr_par.data = [];
    
    load('../models/random_0.07_norm.mat');
    m_parameters(:,:,end) = m_parameters(:,:,end) + 0.8e10 * signal;% .* pattern';
    
    % usr_par.kernel.imfilter.source = fspecial('gaussian', [1 1], 1);
    % usr_par.kernel.imfilter.source = fspecial('gaussian',[75 75], 30);
    % usr_par.kernel.imfilter.source = fspecial('gaussian',[40 40], 20);
    % usr_par.kernel.imfilter.source = fspecial('gaussian',[20 20], 10);
    % usr_par.kernel.imfilter.structure = usr_par.kernel.imfilter.source;    
    usr_par.kernel.sigma.source = [1e-3 1e-3];
    usr_par.kernel.sigma.structure = usr_par.kernel.sigma.source;
    
    [usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);
    
    m_parameters = map_m_to_parameters( map_parameters_to_m(m_parameters, usr_par ) , usr_par );
    m_parameters(:,:,end) = sqrt( m_parameters(:,:,end) ./ rho );
    
%     plot_models( m_parameters, n_basis_fct, array, [0 7 3700 4300] );
%     plot_models( m_parameters, n_basis_fct, array, [0 0 0 0] );

    plot_models_poster( m_parameters, n_basis_fct, array, [0 7 3700 4300] );
    
end



end


