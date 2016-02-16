
clear all





u_0_h = load('~/Desktop/runs/2016_start/data/data_16_ref_0_1h1g_iugg.mat');
u_0_h1g = load('~/Desktop/runs/2016_start/data/data_16_ref_0_1h1g_iugg_newpara.mat');
load('~/Desktop/runs/2016_start/data/array_16_ref.mat');

% u_0_h = load('~/Desktop/runs/2016_start/data/data_16_ref_0_1h_homog_small.mat');
% u_0_h1g = load('~/Desktop/runs/2016_start/data/data_16_ref_0_1h1g_homog_small.mat');
% u_0_h1g = load('~/Desktop/runs/2016_start/data/data_16_ref_0_1h1g_iugg_small.mat');
% load('~/Desktop/runs/2016_start/data/array_16_ref_small.mat');


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
    
%     [misfit( (i-1)*n_rec + 1 : i*n_rec ),~] = make_adjoint_sources_inversion(first(indices,:), second(indices,:), t, veldis, 'log_amplitude_ratio', src, rec);
%     [misfit( (i-1)*n_rec + 1 : i*n_rec ),~] = make_adjoint_sources_inversion(first(indices,:), second(indices,:), t, veldis, 'amplitude_difference', src, rec);
%     [misfit( (i-1)*n_rec + 1 : i*n_rec ),~] = make_adjoint_sources_inversion(first(indices,:), second(indices,:), t, veldis, 'cc_time_shift', src, rec);
    [misfit( (i-1)*n_rec + 1 : i*n_rec ),~] = make_adjoint_sources_inversion(first(indices,:), second(indices,:), t, veldis, 'waveform_difference', src, rec);
    
    
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
% plot_recordings_windows(first(index,:),t,veldis,'k',true,left(index),right(index));
% plot_recordings_windows(second(index,:),t,veldis,'g',true,left(index),right(index));
% plot_recordings(first(index,:),t,veldis,'b',true);
% plot_recordings(second(index,:),t,veldis,'r',true);
