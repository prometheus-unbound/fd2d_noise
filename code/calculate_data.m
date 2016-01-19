
tic

clear all
% close all


% mode = 'local';
mode = 'monch';
% mode = 'euler';
% mode = 'brutus';

use_mex = 'yes';

addpath(genpath('../'))


%% set up model
[Lx,Lz,nx,nz,dt,nt,order,model_type,~,n_basis_fct] = input_parameters();
[X,Z,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);
[width] = absorb_specs();


%% get source and material
[source_dist, spectrum] = make_noise_source();
[mu, rho] = define_material_parameters(nx,nz,model_type);

% mu = mu + 1e9;
% source_dist = source_dist + 1;

if(model_type==666)
    
    A = imread('../models/rand_20_one_sided1.png');
    mu = mu + 5.0e9 * flipud( abs((double(A(:,:,1))-255)/max(max(abs(double(A(:,:,1))-255)))) )';
    mu = mu - 5.0e9 * flipud( abs((double(A(:,:,2))-255)/max(max(abs(double(A(:,:,2))-255)))) )';
    
elseif(model_type==888)
    
    cali = load('california.mat');
    mu = (cali.v).^2 .* rho;
    
end


%% specify output behaviour
output_specs
if(strcmp(mode,'cluster'))
    make_plots = 'no';
end


%% define receiver array
nr_x = 4;
nr_z = 4;
array = zeros(nr_x*nr_z,2);
for i = 1:nr_x
    for j = 1:nr_z        
        array( (i-1)*nr_x + j, 1 ) = 0.9e6 + ( i-1 )*0.25e6;
        array( (i-1)*nr_z + j, 2 ) = 0.6e6 + ( j-1 )*0.25e6;
        
%         array( (i-1)*nr_x + j, 1 ) = 0.7e6 + ( i-1 )*0.2e6;
%         array( (i-1)*nr_z + j, 2 ) = 0.7e6 + ( j-1 )*0.2e6;
        
%         array( (i-1)*nr_x + j, 1 ) = 0.3e6 + ( i-1 )*1.4e6/(nr_x-1);
%         array( (i-1)*nr_z + j, 2 ) = 0.3e6 + ( j-1 )*1.4e6/(nr_z-1);
    end
end


%% small test array, only two receivers close to each other
% array = zeros(2,2);
% array(1,1) = 2.5e4;
% array(2,1) = 3.5e4;
% array(:,2) = 3.0e4;

% array = zeros(2,2);
% array(1,1) = 1.2e5;
% array(2,1) = 2.8e5;
% array(:,2) = 2.0e5;


%% California setup
% array(:,1) = cali.rec_x(1:40);
% array(:,2) = cali.rec_z(1:40);


%% select receivers that will be reference stations
ref_stat = array; %(1,:);
n_ref = size(ref_stat,1);
n_rec = size(array,1)-1;


%% plot configuration
if( strcmp(make_plots,'yes') )
    figure
    hold on
    plot(array(:,1),array(:,2),'o')
    plot(ref_stat(:,1),ref_stat(:,2),'x')
    xlim([0 Lx])
    ylim([0 Lz])
    drawnow
    axis square
    
    return

    plot_model    
end


%% start matlabpool and set up path
if( ~strcmp(mode,'local') )
    addpath(genpath('../'))
    parobj = start_cluster(mode,'', n_ref);
end


%% calculate correlations
t = -(nt-1)*dt:dt:(nt-1)*dt;
c_it = zeros(n_ref,n_rec,length(t));

fprintf('\n')

parfor i = 1:n_ref
    
    if( strcmp(verbose,'yes') )
        fprintf('reference station: %i\n',i)
    end
    
    % each reference station will act as a source once
    src = ref_stat(i,:);
    rec = array( find(~ismember(array,src,'rows') ) , :);
    
    
    fprintf('%i: calculate Green function\n',i)
    if( strcmp(use_mex,'yes') )
        [G_2] = run_forward1_green_mex(mu, rho, src, 1);
    else
        [G_2] = run_forward1_green(mu, rho, src, 1);
    end
    
    fprintf('%i: calculate correlation\n',i)
    if( strcmp(use_mex,'yes') )
        [c_it(i,:,:)] = run_forward2_correlation_mex(mu, rho, G_2, spectrum, source_dist, rec, 1, 0);
    else
        [c_it(i,:,:)] = run_forward2_correlation(mu, rho, G_2, spectrum, source_dist, rec, 1, 0);
    end
    
    fprintf('%i: done\n',i)
    
end


% reorganize correlation vector
c_data = zeros(n_ref*n_rec,length(t));
for i = 1:n_ref
    c_data( (i-1)*n_rec + 1 : i*n_rec, :) = c_it(i,:,:);
end


%% plot data
if( strcmp(make_plots,'yes') )
    figure
    plot_recordings(c_data,t,'vel','k-',true);
    legend('data')
end


%% save array and data for inversion
save( sprintf('../output/interferometry/array_%i_ref.mat',n_ref), 'array', 'ref_stat')
save( sprintf('../output/interferometry/data_%i_ref_%i.mat',n_ref,n_basis_fct), 'c_data', 't')


%% close matlabpool and clean up path
if( ~strcmp(mode,'local') )
    delete(parobj)
    rmpath(genpath('../'))
end


toc
