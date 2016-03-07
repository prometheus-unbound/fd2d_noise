
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
% du = 2 * rand(1,length(u.t)) - 0.5;
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
% usr_par.measurement = 'amplitude_difference';
% usr_par.veldis = 'dis';
% 
% 
% [dcheck, dcheck_struct] = optlib_adjoint_stf_check(u.c_data(index,:),u0.c_data(index,:),du,u.t,src,rec,-10,-2,0.01,usr_par);
% return
% keyboard
% clear all



%% check gradient
[Lx,Lz,nx,nz,dt,nt,order,model_type,source_type,n_basis_fct] = input_parameters();

usr_par.network = load('../output/interferometry/array_2_ref.mat');
usr_par.data = load('../output/interferometry/data_2_ref_0.mat');

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


usr_par.m0 = m;

dm = 0.0 * m_parameters;
dm(:,:,1:end) = dm(:,:,1:end) + 0.1;

[absbound] = init_absbound();

%%% TEST %%%
x_source_r = 3.0e4;
z_source_r = 3.0e4;
radius = 2.0e4;
thickness = 5e3;
[Lx, Lz, nx, nz] = input_parameters();
[X, Z] = define_computational_domain(Lx,Lz,nx,nz);

R = ( (X-x_source_r).^2 + (Z-z_source_r).^2 ).^(1/2);

% mesh( double(R > (radius-thickness/2) & R < (radius+thickness/2) ) )
% return 

% dm(:,:,1) = dm(:,:,1) .* double(R > (radius-thickness/2) & R < (radius+thickness/2) );
%%% TEST END %%%


for i = 1:size(dm,3)
   dm(:,:,i) =  double( absbound == 1 ) .* dm(:,:,i); 
end


[dcheck, dcheck_struct] = optlib_derivative_check( m, reshape(dm,[],1), -10, -1, 1, usr_par);

