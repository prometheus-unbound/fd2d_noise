
function [misfit] = compare_influences( first, second, t, array, ref_stat, plotten )


n_ref = size(ref_stat,1);
n_rec = size(array,1)-1;
distances = zeros(n_ref*n_rec,1);


% veldis = 'vel';
veldis = 'dis';

misfit = zeros(n_rec*n_ref, 4);
for i = 1:n_ref
       
    % each reference station will act as a source once
    src = ref_stat(i,:);
    rec = array( find(~ismember(array,src,'rows') ) , :);
    
    % calculate distance vector
    distances( (i-1)*n_rec + 1 : i*n_rec , 1 ) = sqrt( (src(1,1) - rec(:,1)).^2 + (src(1,2) - rec(:,2)).^2 );
    
    % calculate misfit
    indices = (i-1)*n_rec + 1 : i*n_rec; 
    misfit( (i-1)*n_rec + 1 : i*n_rec, 1 ) = make_adjoint_sources_inversion(first(indices,:), second(indices,:), t, veldis, 'log_amplitude_ratio', src, rec);
    misfit( (i-1)*n_rec + 1 : i*n_rec, 2 ) = make_adjoint_sources_inversion(first(indices,:), second(indices,:), t, veldis, 'amplitude_difference', src, rec);
    misfit( (i-1)*n_rec + 1 : i*n_rec, 3 ) = make_adjoint_sources_inversion(first(indices,:), second(indices,:), t, veldis, 'cc_time_shift', src, rec);
    misfit( (i-1)*n_rec + 1 : i*n_rec, 4 ) = make_adjoint_sources_inversion(first(indices,:), second(indices,:), t, veldis, 'waveform_difference', src, rec);
    
    
end


left = distances/4000.0 - 27.0;
right = distances/4000.0 + 27.0;

i_zero = find( t==0 );
i_left_neg = find( left<0 );
left(i_left_neg) = t(i_zero+1);

if( right > t(end) )
    right = t(end);
end


fprintf('loga:  %8.4f,   max: %8.4f,   mean: %8.4f\n', sum(abs(misfit(:,1))), max(abs(misfit(:,1))), mean(abs(misfit(:,1))) )
fprintf('amp:   %8.4f,   max: %8.4f,   mean: %8.4f\n', sum(abs(misfit(:,2))), max(abs(misfit(:,2))), mean(abs(misfit(:,2))) )
fprintf('cc:    %8.4f,   max: %8.4f,   mean: %8.4f\n', sum(abs(misfit(:,3))), max(abs(misfit(:,3))), mean(abs(misfit(:,3))) )
fprintf('wd:    %8.4f,   max: %8.4f,   mean: %8.4f\n', sum(abs(misfit(:,4))), max(abs(misfit(:,4))), mean(abs(misfit(:,4))) )

if( strcmp(plotten,'yes') )
    figure
    % % index = 165:180;
    index = 1:n_ref*n_rec;
    % plot_recordings_windows(first(index,:),t,veldis,'k',true,left(index),right(index));
    % plot_recordings_windows(second(index,:),t,veldis,'g',true,left(index),right(index));
    plot_recordings(first(index,:),t,veldis,'b',true);
    plot_recordings(second(index,:),t,veldis,'r',true);
end


end