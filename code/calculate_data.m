
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use_mex = 'yes';
% 'yes' (if startup.m indicates a successful compilation)
% 'no' (default)

% small test array for gradient validation
% array = zeros(2,2);
% array(1,1) = 2.5e4;
% array(2,1) = 3.5e4;
% array(:,2) = 3.0e4;

% array for kernel computation
array = zeros(2,2);
array(1,1) = 1.25e5;
array(2,1) = 2.75e5;
array(:,2) = 2.0e5;

% select receivers that will be reference stations
ref_stat = array(1,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get configuration
[~,~,nx,nz,dt,nt,~,model_type,source_type,~,make_plots] = input_parameters();


% get source and material
noise_source = make_noise_source('no');
structure = define_material_parameters('no');


% for gradient test
% noise_source.distribution = noise_source.distribution + rand(nx,nz);
% structure.mu = structure.mu + 1e9;


% plot model with array configuration
if( strcmp(make_plots,'yes') )
    plot_models( sqrt(structure.mu./structure.rho), noise_source.distribution, array, [0 0 0 0]);
end


if( strcmp( use_mex, 'no' ) )
    ! rm ../code_mex_functions/run*
    ! cp ../code/run1_forward_green.m ../code/mex_functions/run1_forward_green_mex.m
    ! cp ../code/run2_forward_correlation.m ../code/mex_functions/run2_forward_correlation_mex.m
    ! cp ../code/run3_adjoint.m ../code/mex_functions/run3_adjoint_mex.m
end


% loop over reference stations
n_ref = size(ref_stat,1);
n_rec = size(array,1)-1;
t = -(nt-1)*dt:dt:(nt-1)*dt;
nt = length(t);
c_iref = zeros(n_ref,n_rec,length(t));
fprintf('\n')

tic
for i_ref = 1:n_ref
        
    src = ref_stat(i_ref,:);
    rec = array( find(~ismember(array,src,'rows') ) , :);
    
    fprintf( 'ref %i: calculate Green function\n', i_ref )    
    % if( ~exist(['../output/G_fft_ref_' num2str(i_ref) '_model_' num2str(model_type) '.mat'], 'file') )
        G_fft = run1_forward_green_mex( structure, src, 0 );
    %     parsave( ['../output/G_fft_ref_' num2str(i_ref) '_model_' num2str(model_type) '.mat'], G_fft )
    % else
    %     G_fft = parload( ['../output/G_fft_ref_' num2str(i_ref) '_model_' num2str(model_type) '.mat'] );
    % end
    
    fprintf( 'ref %i: calculate correlation\n', i_ref )    
    c_iref(i_ref,:,:) = run2_forward_correlation_mex( structure, noise_source, G_fft, rec, 0 );
    
    fprintf( 'ref %i: done\n', i_ref )
    
end
toc


% reorganize correlation vector (for parfor-users)
c_data = zeros(n_ref*n_rec,length(t));
for i_ref = 1:n_ref
    c_data( (i_ref-1)*n_rec + 1 : i_ref*n_rec, :) = c_iref(i_ref,:,:);
end


% plot data
if( strcmp(make_plots,'yes') )
    figure
    plot_recordings(c_data, t, 'vel', 'k', true);
    legend('data')
end


% try path
current_path = pwd;
path_index = strfind(current_path, 'fd2d_noise');


% save array and data for inversion
save( [ current_path(1:path_index+9) '/output/array_' num2str(n_ref) '_ref.mat' ], 'array', 'ref_stat')
save( [ current_path(1:path_index+9) '../output/data_' num2str(n_ref) '_ref_model_' model_type '_source_' source_type '.mat' ], 'c_data', 't')

