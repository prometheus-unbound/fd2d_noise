
clear all


% folder_1 = '~/Desktop/runs/inversion/data/';
% 
% u_h = load([folder_1 'data_16_ref_uniform_homog_structure.mat']);
% u_1 = load([folder_1 'data_16_ref_uniform_structure_1.mat']);
% u_2 = load([folder_1 'data_16_ref_uniform_structure_2.mat']);
% u_3 = load([folder_1 'data_16_ref_uniform_structure_3.mat']);
% 
% b3_1 = load([folder_1 'data_16_ref_blob3_structure_1.mat']);
% b3_h = load([folder_1 'data_16_ref_blob3_homog_structure.mat']);
% 
% ub3_1 = load([folder_1 'data_16_ref_uniform_blob3_structure_1.mat']);
% ub3_h = load([folder_1 'data_16_ref_uniform_blob3_homog_structure.mat']);
% ub3_r5 = load([folder_1 'data_16_ref_uniform_blob3_rand_5.mat']);
% ub3_r10 = load([folder_1 'data_16_ref_uniform_blob3_rand_10.mat']);
% ub3_r10m = load([folder_1 'data_16_ref_uniform_blob3_rand_10_muchos.mat']);
% ub3_r10d = load([folder_1 'data_16_ref_uniform_blob3_rand_10_demasiados.mat']);
% ub3_r20 = load([folder_1 'data_16_ref_uniform_blob3_rand_20.mat']);
% ub3_r30 = load([folder_1 'data_16_ref_uniform_blob3_rand_30.mat']);
% ub3_r40 = load([folder_1 'data_16_ref_uniform_blob3_rand_40.mat']);
% 
% u2b_3_1_1 = load([folder_1 'data_16_ref_uniform_2blob_3_1_structure_1.mat']);
% 
% ub20_1 = load([folder_1 'data_16_ref_uniform_blob20_structure_1.mat']);
% ub20_h = load([folder_1 'data_16_ref_uniform_blob20_homog_structure.mat']);
% 
% ub100_1 = load([folder_1 'data_16_ref_uniform_blob100_structure_1.mat']);
% ub100_2 = load([folder_1 'data_16_ref_uniform_blob100_structure_2.mat']);
% ub100_3 = load([folder_1 'data_16_ref_uniform_blob100_structure_3.mat']);
% ub100_h = load([folder_1 'data_16_ref_uniform_blob100_homog_structure.mat']);
% 
% ub1000_1 = load([folder_1 'data_16_ref_uniform_blob1000_structure_1.mat']);
% ub1000_h = load([folder_1 'data_16_ref_uniform_blob1000_homog_structure']);



% folder_1 = '~/Desktop/runs/inversion/data/';
% 
% u_h_0 = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_0_uniform_homogeneous.mat');
% u_h_5 = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_5_uniform_homogeneous.mat');
% u_h_10 = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_10_uniform_homogeneous.mat');
% u_h_20 = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_20_uniform_homogeneous.mat');
% 
% u1_h_0 = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_0_uniform1_homogeneous.mat');
% u2_h_0 = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_0_uniform2_homogeneous.mat');
% 
% u1b_h_0 = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_0_uniform_1gaussian_homogeneous.mat');
% u2b_h_0 = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_0_uniform_2gaussian_homogeneous.mat');
% 
% load([folder_1 'array_16_ref.mat'])
% t = u_h_0.t;



folder_1 = '~/Desktop/runs/inversion_newest/data/';

% u_0_h = load('~/Desktop/runs/inversion_newest/data/data_16_ref_0_1h.mat');
% u_0_h1g = load('~/Desktop/runs/inversion_newest/data/data_16_ref_0_1h1g.mat');
% u_0_h2g = load('~/Desktop/runs/inversion_newest/data/data_16_ref_0_2h2g.mat');
% u_0_h2g_n = load('~/Desktop/runs/inversion_newest/data/data_16_ref_0_2h2g_lr_nover.mat');
% u_0_r5os = load('~/Desktop/runs/inversion_newest/data/data_16_ref_0_rand_5_one_sided1.mat');
% u_0_r10os = load('~/Desktop/runs/inversion_newest/data/data_16_ref_0_rand_10_one_sided1.mat');
% u_0_r15os = load('~/Desktop/runs/inversion_newest/data/data_16_ref_0_rand_15_one_sided1.mat');
% u_0_r20os = load('~/Desktop/runs/inversion_newest/data/data_16_ref_0_rand_20_one_sided1.mat');
% 
% u_0_h1g_p_fs1 = load('~/Desktop/runs/inversion_newest/data/point/data_16_ref_0_1h1g_point1_freqsamp1.mat');
% u_0_h1g_p_fs5 = load('~/Desktop/runs/inversion_newest/data/point/data_16_ref_0_1h1g_point1_freqsamp5.mat');
% u_0_h_p_fs1 = load('~/Desktop/runs/inversion_newest/data/point/data_16_ref_0_1h_point1_freqsamp1.mat');
% u_0_h_p_fs5 = load('~/Desktop/runs/inversion_newest/data/point/data_16_ref_0_1h_point1_freqsamp5.mat');

u_0_h = load('~/Desktop/runs/2016_start/data/data_16_ref_0_1h_homog_small.mat');
u_0_h1g = load('~/Desktop/runs/2016_start/data/data_16_ref_0_1h1g_iugg_small.mat');

load('~/Desktop/runs/2016_start/data/array_16_ref_small.mat');

% load([folder_1 'array_16_ref.mat'])
t = u_0_h.t;



n_ref = size(ref_stat,1);
n_rec = size(array,1)-1;
distances = zeros(n_ref*n_rec,1);
first = zeros(n_ref*n_rec,length(t));
second = zeros(n_ref*n_rec,length(t));


% veldis = 'vel';
veldis = 'dis';

f_min = 1/15 - 0.01;
f_max = 1/15 + 0.01;

misfit = 0;
for i = 1:n_ref
       
    % each reference station will act as a source once
    src = ref_stat(i,:);
    rec = array( find(~ismember(array,src,'rows') ) , :);
    
    % calculate distance vector
    distances( (i-1)*n_rec + 1 : i*n_rec , 1 ) = sqrt( (src(1,1) - rec(:,1)).^2 + (src(1,2) - rec(:,2)).^2 );
    
    % calculate misfit
    indices = (i-1)*n_rec + 1 : i*n_rec;
    
    % first(indices,:) = u2_h_0.c_data( indices , : );
    % second(indices,:) = u2b_h_0.c_data( indices , : );
    
    first(indices,:) = u_0_h.c_data( indices , : );
    second(indices,:) = u_0_h1g.c_data( indices , : );
    
%     first_save = first;
%     second_save = second;  
%     first(indices,:) = filter_correlations( first(indices,:), t, f_min, f_max );
%     second(indices,:) = filter_correlations( second(indices,:), t, f_min, f_max );
    
    [misfit( (i-1)*n_rec + 1 : i*n_rec ),~] = make_adjoint_sources_inversion(first(indices,:), second(indices,:), t, veldis, 'log_amplitude_ratio', src, rec);
%     [misfit( (i-1)*n_rec + 1 : i*n_rec ),~] = make_adjoint_sources_inversion(first(indices,:), second(indices,:), t, veldis, 'amplitude_difference', src, rec);
%     [misfit( (i-1)*n_rec + 1 : i*n_rec ),~] = make_adjoint_sources_inversion(first(indices,:), second(indices,:), t, veldis, 'cc_time_shift', src, rec);
%     [misfit( (i-1)*n_rec + 1 : i*n_rec ),~] = make_adjoint_sources_inversion(first(indices,:), second(indices,:), t, veldis, 'waveform_difference', src, rec);
    
    
end


left = distances/4000.0 - 27.0;
right = distances/4000.0 + 27.0;

i_zero = find( t==0 );
i_left_neg = find( left<0 );
left(i_left_neg) = t(i_zero+1);

if( right > t(end) )
    right = t(end);
end


fprintf('%f\n',sum(abs(misfit)))
% fprintf('%f',max(abs(misfit)))


% % index = 165:180;
index = 1:n_ref*n_rec;
plot_recordings_windows(first(index,:),t,veldis,'k',true,left(index),right(index));
plot_recordings_windows(second(index,:),t,veldis,'g',true,left(index),right(index));
% plot_recordings(first(index,:),t,veldis,'b',true);
% plot_recordings(second(index,:),t,veldis,'r',true);
