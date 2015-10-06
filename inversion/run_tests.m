
clear all
close all


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
% usr_par.measurement = 'waveform_difference';
% usr_par.veldis = 'dis';
% 
% 
% [dcheck, dcheck_struct] = optlib_adjoint_stf_check(u.c_data(index,:),u0.c_data(index,:),du,u.t,src,rec,-10,-2,0.01,usr_par);
% keyboard
% clear all



%% check adjoint state
% [~,~,~,~,~,nt] = input_parameters();
% m = make_noise_source();
% df = m;
% 
% df_time_dependent = 0.0 * repmat( df, 1, 1, 2*nt-1 );
% for n = 1:size(df_time_dependent,3)
% 
%     if( mod(n,5) == 2 )
%         df_time_dependent(:,:,n) = 1.0;
%     else
%         df_time_dependent(:,:,n) = 0.0;
%     end
%     
% end
% 
% usr_par = usr_par_init_default_parameters_lbfgs([]);
% usr_par.debug.switch = 'yes';
% 
% [dcheck, dcheck_struct] = optlib_adjoint_state_check( reshape(m,[],1), df_time_dependent, -10, -2, 1, usr_par);
% keyboard
% clear all



%% check gradient
% m = reshape( make_noise_source(), [], 1);
% dm = 0.1 * m;

% [Lx,Lz,nx,nz,dt,nt,order,model_type,source_type,n_basis_fct] = input_parameters();
% [mu,~] = define_material_parameters(nx,nz,model_type);
% [usr_par] = usr_par_init_default_parameters_lbfgs([]);
% m = map_parameters_to_m(mu,usr_par);
% dm = m + 0.1;

% [dcheck, dcheck_struct] = optlib_derivative_check( m, dm, -3, 1, 0.5, usr_par);

