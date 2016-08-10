
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use_mex = 'no';
% 'yes' (if startup.m indicates a successful compilation)
% 'no' (default)

% small test array for gradient validation
% array = zeros(2,2);
% array(1,1) = 2.5e4;
% array(2,1) = 3.5e4;
% array(:,2) = 3.0e4;

% array for kernel computation
array = zeros(2,2);
array(1,1) = 1.3e5;
array(2,1) = 2.7e5;
array(:,2) = 2.0e5;

% select receivers that will be reference stations
ref_stat = array(1,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% check path and mex files
fd2d_path = fd2d_path();
check_mex_files( use_mex );    


% get configuration
[~,~,nx,nz,dt,nt,~,model_type,source_type,~,make_plots] = input_parameters();


% get source and material
noise_source = make_noise_source('no');
structure = define_material_parameters('no');


% for gradient test
% noise_source.distribution = noise_source.distribution + rand(nx,nz);
% structure.mu = structure.mu + 1e9;


% plot model with array configuration
% if( strcmp(make_plots,'yes') )
%     plot_models( sqrt(structure.mu./structure.rho), noise_source.distribution, array, [0 7 0 0]);
% end


% loop over reference stations
n_ref = size(ref_stat,1);
n_rec = size(array,1)-1;
t = -(nt-1)*dt:dt:(nt-1)*dt;
nt = length(t);
correlations = zeros(n_ref,n_rec,length(t));
fprintf('\n')

tic
for i_ref = 1:n_ref
        
    src = ref_stat(i_ref,:);
    rec = array( find(~ismember(array,src,'rows') ) , :);
    
    fprintf( 'ref %i: calculate Green function\n', i_ref )
    if( ~exist( filename('G_fft', i_ref), 'file' ) )
        G_fft = run1_forward_green_mex( structure, src, 0 );
        parsave( filename('G_fft', i_ref), G_fft, [] )
    else
        G_fft = parload( filename('G_fft', i_ref) );
    end    
    
    fprintf( 'ref %i: calculate correlation\n', i_ref )    
    [correlations(i_ref,:,:), C] = run2_forward_correlation_mex( structure, noise_source, G_fft, src, rec, 1 );
    
    fprintf( 'ref %i: done\n', i_ref )
    
end
toc


% plot data
if( strcmp(make_plots,'yes') )
    figure
    plot_recordings(correlations, t, 'k', true);
    legend('data')
end


% save array and data for inversion
save( filename('array', n_ref), 'array', 'ref_stat')
save( filename('correlations', n_ref), 'correlations')

