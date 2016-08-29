
clear all
% close all


%% check adjoint source time functions
% folder_1 = '~/Desktop/runs/inversion/data/';
% 
% load([folder_1 'array_16_ref.mat'])
% u = load([folder_1 'data_16_ref_uniform_homog_structure.mat']);
% u0 = load([folder_1 'data_16_ref_uniform_structure_1.mat']);
% 
% % choose du randomly
% du = 2 * (rand(1,length(u.t)) - 0.5);
% 
% % choose one reference station and one receiver
% i_ref = 2;
% rec_id = 11;
% 
% src = ref_stat(i_ref,:);
% rec = array( find(~ismember(array,src,'rows') ) , :);
% index = (i_ref-1)*size(rec,1) + 1 : i_ref*size(rec,1);
% index = index(rec_id);
% rec = rec(rec_id,:);
% 
% % specify type of measurement and units of correlation
% % usr_par.measurement = 'waveform_difference';
% % usr_par.measurement = 'log_amplitude_ratio';
% % usr_par.measurement = 'amplitude_difference';
% usr_par.veldis = 'dis';
% 
% 
% [dcheck, dcheck_struct] = optlib_check_adjoint_stf( u.c_data(index,:), u0.c_data(index,:), du, u.t, src, rec, -10, -2, 0.01, usr_par );
% 
% 
% keyboard
% clear all




%% check gradient
[Lx,Lz,nx,nz,dt,nt,order,model_type,source_type,n_basis_fct] = input_parameters();

usr_par.network = load('../output/interferometry/array_1_ref.mat');
usr_par.data = load('../output/interferometry/data_1_ref_0.mat');
usr_par.type = 'joint';
usr_par.kernel.weighting = 0.5;

usr_par.measurement.source = 'waveform_difference';
usr_par.measurement.structure = 'waveform_difference';
        
usr_par.kernel.sigma.source = [1e-3 1e-3];
usr_par.kernel.sigma.structure = usr_par.kernel.sigma.source;

[usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);


% set up initial model
if( n_basis_fct == 0)
    m_parameters = zeros(nx, nz, 2);
    m_parameters(:,:,1) = make_noise_source(source_type, n_basis_fct);
else
    m_parameters = zeros(nx, nz, n_basis_fct+1);
    m_parameters(:,:,1:n_basis_fct) = make_noise_source(source_type, n_basis_fct);
end
m_parameters(:,:,end) = define_material_parameters(nx,nz,model_type);


% convert to optimization variable
m = map_parameters_to_m(m_parameters,usr_par);


% initial model for regularization
usr_par.m0 = m;


% set up test vector
dm = rand( numel(m), 1 );


%%% TEST %%%
% x_source_r = 3.0e4;
% z_source_r = 3.0e4;
% radius = 2.0e4;
% thickness = 5e3;
% [Lx, Lz, nx, nz] = input_parameters();
% [X, Z] = define_computational_domain(Lx,Lz,nx,nz);
% 
% R = ( (X-x_source_r).^2 + (Z-z_source_r).^2 ).^(1/2);
% 
% mesh( double(R > (radius-thickness/2) & R < (radius+thickness/2) ) )
% return 
% 
% dm(:,:,1) = dm(:,:,1) .* double(R > (radius-thickness/2) & R < (radius+thickness/2) );
%%% TEST END %%%


% make sure that test vector is consistent with setup
if( strcmp( usr_par.type, 'source' ) && ~isempty(find( dm(end-nx*nz+1:end, 1), 1 )) )
    fprintf('\nchanged dm accordingly!\n')
    dm( end - nx*nz + 1 : end, 1 ) = 0 * dm( end - nx*nz + 1 : end, 1 );
elseif( strcmp( usr_par.type, 'structure' ) && ~isempty(find( dm(1:end-nx*nz, 1), 1 )) )
    fprintf('\nchanged dm accordingly!\n')
    dm( 1:end-nx*nz, 1 ) = 0 * dm( 1:end-nx*nz, 1 );
end


[absbound] = reshape( init_absbound() , [], 1 );
for i = 1:size(dm)/(nx*nz)
   dm( ((i-1)*nx*nz+1):i*nx*nz, 1 ) = double( absbound == 1 ) .* dm( ((i-1)*nx*nz+1):i*nx*nz, 1 ); 
end


[dcheck, dcheck_struct] = optlib_check_derivative( m, reshape(dm,[],1), -10, 0, 1, usr_par );

return
keyboard
clear all




%% check adjoint source time functions for hessian vector products
% folder_1 = '~/Desktop/runs/inversion/data/';
% 
% load([folder_1 'array_16_ref.mat'])
% u = load([folder_1 'data_16_ref_uniform_homog_structure.mat']);
% u0 = load([folder_1 'data_16_ref_uniform_structure_1.mat']);
% 
% % choose du randomly
% du1 = 1e-3 * 2 * (rand(1,length(u.t)) - 0.5);
% du2 = 1e-3 * 2 * (rand(1,length(u.t)) - 0.5);
% 
% % choose one reference station and one receiver
% i_ref = 2;
% rec_id = 11;
% 
% src = ref_stat(i_ref,:);
% rec = array( find(~ismember(array,src,'rows') ) , :);
% index = (i_ref-1)*size(rec,1) + 1 : i_ref*size(rec,1);
% index = index(rec_id);
% rec = rec(rec_id,:);
% 
% % specify type of measurement and units of correlation
% % usr_par.measurement = 'waveform_difference';
% % usr_par.measurement = 'log_amplitude_ratio';
% usr_par.measurement = 'amplitude_difference';
% usr_par.veldis = 'dis';
% 
% 
% [dcheck, dcheck_struct] = optlib_check_adjoint_stf_hessian( u.c_data(index,:), u0.c_data(index,:), du1, du2, u.t, src, rec, -10, -1, 1, usr_par );
% 
% 
% keyboard
% clear all




%% check du - forward and adjoint
% [Lx,Lz,nx,nz,dt,nt,order,model_type,source_type,n_basis_fct] = input_parameters();
% [X,Z] = define_computational_domain(Lx,Lz,nx,nz);
% 
% usr_par.network = load('../output/interferometry/array_1_ref.mat');
% usr_par.data = load('../output/interferometry/data_1_ref_0.mat');
% [usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);
% 
% 
% % set up initial model
% if( n_basis_fct == 0)
%     m_parameters = zeros( nx, nz, 2 );
%     m_parameters(:,:,1) = make_noise_source( source_type, n_basis_fct );
% else
%     m_parameters = zeros( nx, nz, n_basis_fct+1 );
%     m_parameters(:,:,1:n_basis_fct) = make_noise_source( source_type, n_basis_fct );
% end
% m_parameters(:,:,end) = define_material_parameters( nx, nz, model_type );
% 
% 
% % convert to optimization variable
% m = map_parameters_to_m( m_parameters, usr_par );
% 
% 
% % initial model for regularization
% usr_par.m0 = m;
% 
% 
% % set up test vectors
% dm = 0 * m;
% dm( nx*nz+1:end, 1 ) = rand( nx*nz, 1 );
% 
% 
% % set test vectors in absorbing boundary region to zero
% [absbound] = reshape( init_absbound() , [], 1 );
% for i = 1:size(dm)/(nx*nz)
%    dm( ((i-1)*nx*nz+1):i*nx*nz, 1 ) = double( absbound == 1 ) .* dm( ((i-1)*nx*nz+1):i*nx*nz, 1 ); 
% end
% 
% 
% [dcheck, dcheck_struct] = optlib_check_du( m, dm, -10, -1, 1, usr_par );
% 
% return
% keyboard
% clear all




%% check Hessian vector product
[Lx,Lz,nx,nz,dt,nt,order,model_type,source_type,n_basis_fct] = input_parameters();
[X,Z] = define_computational_domain(Lx,Lz,nx,nz);

usr_par.network = load('../output/interferometry/array_1_ref_testing_small.mat');
usr_par.data = load('../output/interferometry/data_1_ref_0_testing_small.mat');

% usr_par.network = load('../output/interferometry/array_1_ref.mat');
% usr_par.data = load('../output/interferometry/data_1_ref_0.mat');

usr_par.use_mex = 'no';
usr_par.type = 'source';

usr_par.kernel.weighting = 0.5;
usr_par.measurement.source = 'waveform_difference';
% usr_par.measurement.source = 'log_amplitude_ratio';
% usr_par.measurement.source = 'amplitude_difference';
usr_par.measurement.structure = 'waveform_difference';
% usr_par.measurement.structure = 'log_amplitude_ratio';

usr_par.regularization.alpha = 0;
usr_par.regularization.beta = 0;

[usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);


% set up initial model
if( n_basis_fct == 0)
    m_parameters = zeros( nx, nz, 2 );
    m_parameters(:,:,1) = make_noise_source( source_type, n_basis_fct );
else
    m_parameters = zeros( nx, nz, n_basis_fct+1 );
    m_parameters(:,:,1:n_basis_fct) = make_noise_source( source_type, n_basis_fct );
end
m_parameters(:,:,end) = define_material_parameters( nx, nz, model_type );


% convert to optimization variable
m = map_parameters_to_m( m_parameters, usr_par );


% initial model for regularization
usr_par.m0 = m;


% set up test vectors
dm1 = rand( numel(m), 1 );
dm2 = rand( numel(m), 1 );


% make sure that test vectors are consistent with setup
if( strcmp( usr_par.type, 'source' ) && ~isempty(find( dm1(end-nx*nz+1:end, 1), 1 )) )
    fprintf('\nchanged dm1 accordingly!\n')
    dm1( end - nx*nz + 1 : end, 1 ) = 0 * dm1( end - nx*nz + 1 : end, 1 );
elseif( strcmp( usr_par.type, 'structure' ) && ~isempty(find( dm1(1:end-nx*nz, 1), 1 )) )
    fprintf('\nchanged dm1 accordingly!\n')
    dm1( 1:end-nx*nz, 1 ) = 0 * dm1( 1:end-nx*nz, 1 );
end


% set test vectors in absorbing boundary region to zero
[absbound] = reshape( init_absbound() , [], 1 );
for i = 1:size(dm1)/(nx*nz)
   dm1( ((i-1)*nx*nz+1):i*nx*nz, 1 ) = double( absbound == 1 ) .* dm1( ((i-1)*nx*nz+1):i*nx*nz, 1 ); 
   dm2( ((i-1)*nx*nz+1):i*nx*nz, 1 ) = double( absbound == 1 ) .* dm2( ((i-1)*nx*nz+1):i*nx*nz, 1 );
end


tic
[dcheck, dcheck_struct] = optlib_check_hessian( m, dm1, dm2, -9, -1, 1, usr_par );
toc

