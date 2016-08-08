
clear all

addpath(genpath('../'))
[Lx, Lz, nx, nz, dt, nt, ~, model_type, source_type, ~, make_plots] = input_parameters();
[X,Z,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


% get source and material
noise_source = make_noise_source('no');
structure = define_material_parameters('no');


% plot array configuration
if( strcmp(make_plots,'yes') )
    plot_models( sqrt(structure.mu./structure.rho), noise_source.distribution, array, [0 0 0 0]);
end


% loop over reference stations
n_ref = size(ref_stat,1);
n_rec = size(array,1)-1;
t = -(nt-1)*dt:dt:(nt-1)*dt;
c_it = zeros(n_ref,n_rec,length(t));
fprintf('\n')

tic
for i = 1:n_ref
        
    src = ref_stat(i,:);
    rec = array( find(~ismember(array,src,'rows') ) , :);
    
    fprintf( 'ref %i: calculate Green function\n', i )    
    if( ~exist(['../output/G_fft_ref_' num2str(i) '_model_' num2str(model_type) '.mat'], 'file') )
        G_fft = run_forward1_green( structure, src, 0 );
        parsave( ['../output/G_fft_ref_' num2str(i) '_model_' num2str(model_type) '.mat'], G_fft )
    else
        G_fft = parload( ['../output/G_fft_ref_' num2str(i) '_model_' num2str(model_type) '.mat'] );
    end
    
    fprintf( 'ref %i: calculate correlation\n', i )    
    c_it(i,:,:) = run_forward2_correlation( structure, noise_source, G_fft, rec, 0 );    
    
    fprintf( 'ref %i: done\n', i )
    
end
toc


% reorganize correlation vector (for parfor-users)
c_data = zeros(n_ref*n_rec,length(t));
for i = 1:n_ref
    c_data( (i-1)*n_rec + 1 : i*n_rec, :) = c_it(i,:,:);
end


% plot data
if( strcmp(make_plots,'yes') )
    figure
    plot_recordings(c_data, t, 'vel', 'k', true);
    legend('data')
end


% save array and data for inversion
save( sprintf('../output/array_%i_ref.mat',n_ref), 'array', 'ref_stat')
save( sprintf('../output/data_%i_ref_model_%i_source_%s.mat',n_ref,model_type,source_type), 'c_data', 't')

