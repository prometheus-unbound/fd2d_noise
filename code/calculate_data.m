
clear all

addpath(genpath('../'))
[Lx,Lz,nx,nz,dt,nt,~,model_type,source_type,n_basis_fct] = input_parameters();
[X,Z,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

usr_par.cluster = 'monch';
% 'local';
% 'monch';
% 'euler';
% 'brutus';


usr_par.use_mex = 'yes';
% 'yes'
% 'no'


% define reference models for perturbations and mu_0 for structure
usr_par.initial.ref_source = 0;
usr_par.initial.ref_structure = 1;
usr_par.initial.mu_0 = 4.8e10;

% usr_par.kernel.imfilter.source = fspecial('gaussian',[75 75], 30);
% usr_par.kernel.imfilter.source = fspecial('gaussian',[40 40], 20);
% usr_par.kernel.imfilter.source = fspecial('gaussian',[20 20], 10);
% usr_par.kernel.imfilter.source = fspecial('gaussian',[1 1], 1);
% usr_par.kernel.imfilter.structure = usr_par.kernel.imfilter.source;

usr_par.kernel.sigma.source = [5e4 5e4];
usr_par.kernel.sigma.structure = usr_par.kernel.sigma.source;


% define receiver array
nr_x = 4;
nr_z = 4;
array = zeros(nr_x*nr_z,2);
for i = 1:nr_x
    for j = 1:nr_z        
        
        % big - normal config
        array( (i-1)*nr_z + j, 1 ) = 0.9e6 + ( i-1 ) * 0.25e6;
        array( (i-1)*nr_z + j, 2 ) = 0.6e6 + ( j-1 ) * 0.25e6;

        
        % small - normal config
        % array( (i-1)*nr_z + j, 1 ) = 1.8e5 + ( i-1 ) * 0.5e5;
        % array( (i-1)*nr_z + j, 2 ) = 1.2e5 + ( j-1 ) * 0.5e5;

        % small - full coverage
        % array( (i-1)*nr_z + j, 1 ) = 0.8e5 + ( i-1 ) * 0.8e5;
        % array( (i-1)*nr_z + j, 2 ) = 0.8e5 + ( j-1 ) * 0.8e5;

        % small - japan
        % array( (i-1)*nr_z + j, 1 ) = 0.8e5 + ( i-1 ) * 0.8e5;
        % array( (i-1)*nr_z + j, 2 ) = 0.8e5 + ( j-1 ) * 0.4e5;
        
        % small - japan
        % array( nr_x*nr_z + (i-1)*nr_z + j, 1 ) = 0.8e5 + ( i-1 ) * 0.8e5;
        % array( nr_x*nr_z + (i-1)*nr_z + j, 2 ) = 2.8e5 + ( j-1 ) * 0.4e5;

                
        % big - line setup
        %array( (i-1)*nr_z + j, 1 ) = 0.5e6 + ( i-1 ) * 0.1e6;
        % array( (i-1)*nr_z + j, 2 ) = 1.0e6;
        
        % small - line setup
        % array( (i-1)*nr_z + j, 1 ) = 1.0e5 + ( i-1 ) * 0.1e5;
        % array( (i-1)*nr_z + j, 2 ) = 2.0e5;
        
    end
end


% small test array, only two receivers close to each other
% array = zeros(2,2);
% array(1,1) = 2.5e4;
% array(2,1) = 3.5e4;
% array(:,2) = 3.0e4;

% array = zeros(2,2);
% array(1,1) = 1.25e5;
% array(2,1) = 2.75e5;
% array(:,2) = 2.0e5;


% select receivers that will be reference stations
ref_stat = array; %(1,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set necessary fields that might not have been set
usr_par.network = [];
usr_par.data = [];
[usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);


% get source and material
if( n_basis_fct == 0 )
    m_parameters = zeros( nx, nz, 2 );
else
    m_parameters = zeros( nx, nz, n_basis_fct+1 );
end

m_parameters(:,:,1:end-1) = make_noise_source( source_type, n_basis_fct );
m_parameters(:,:,end) = define_material_parameters( nx, nz, model_type );


% m_parameters(:,:,end) = m_parameters(:,:,end) + 1e9;
% if( n_basis_fct == 0 )
%     m_parameters(:,:,1:end-1) = m_parameters(:,:,1:end-1) + rand(nx,nz);
% else
%     m_parameters(:,:,1:end-1) = m_parameters(:,:,1:end-1) + rand(nx,nz,n_basis_fct);
% end


if(model_type==666)
    
    A = imread('../models/random_10_in_lt_array.png');
    m_parameters(:,:,end) = m_parameters(:,:,end) + 5.0e9 * flipud( abs((double(A(:,:,1))-255)/max(max(abs(double(A(:,:,1))-255)))) )';
    m_parameters(:,:,end) = m_parameters(:,:,end) - 5.0e9 * flipud( abs((double(A(:,:,2))-255)/max(max(abs(double(A(:,:,2))-255)))) )';
        
elseif(model_type==888)   
    
    load('../models/random_0.1_norm.mat');
    m_parameters(:,:,end) = m_parameters(:,:,end) + 150.0e9 * signal;
    
end


% convert to m and then back to parameters again, necessary since the smoothing operator is part of the parameterization
m_parameters = map_m_to_parameters( map_parameters_to_m(m_parameters, usr_par ) , usr_par );


% redirect parameters 
source_dist = m_parameters(:,:,1:end-1);
[~, spectrum] = make_noise_source( source_type, n_basis_fct );             % important for n_basis_fct == 0
mu = m_parameters(:,:,end);
[~,rho] = define_material_parameters( nx, nz, model_type );


% specify output behaviour
output_specs
if(strcmp(usr_par.cluster,'cluster'))
    make_plots = 'no';
end


% plot configuration
% if( strcmp(make_plots,'yes') )
%     
%     figure
%     hold on
%     plot( array(:,1), array(:,2), 'o' )
%     plot( ref_stat(:,1), ref_stat(:,2), 'x' )
%     
%     [absorb_width] = absorb_specs();
%     level = [1, 1];
%     plot3( [absorb_width, Lx-absorb_width], [absorb_width, absorb_width], level, 'k--' )
%     plot3( [absorb_width, Lx-absorb_width], [Lz-absorb_width, Lz-absorb_width], level, 'k--' )
%     plot3( [absorb_width, absorb_width], [absorb_width, Lz-absorb_width], level, 'k--' )
%     plot3( [Lx-absorb_width, Lx-absorb_width], [absorb_width, Lz-absorb_width], level, 'k--' )
%     
%     axis image
%     xlim([0 Lx])
%     ylim([0 Lz])
%     drawnow
%     
%     return
%     
%     figure
%     hold on
%     mesh( X, Z, sqrt(mu./rho)' )
%     plot3( array(:,1), array(:,2), 5000 + 0*array(:,2), 'o' )
%     xlim([0 Lx])
%     ylim([0 Lz])
%     axis square
%     colorbar
%     grid on
%     
%     return
% end


% calculate correlations
n_ref = size(ref_stat,1);
n_rec = size(array,1)-1;
t = -(nt-1)*dt:dt:(nt-1)*dt;
c_it = zeros(n_ref,n_rec,length(t));
fprintf('\n')


% start matlabpool
if( ~strcmp(usr_par.cluster,'local') )
    parobj = start_cluster( usr_par.cluster, '', n_ref );
end


tic
parfor i = 1:n_ref
        
    src = ref_stat(i,:);
    rec = array( find(~ismember(array,src,'rows') ) , :);
    
    fprintf( '%i: calculate Green function\n', i )
    G = run_forward1_green_mex( mu, rho, src, rec, 0, [], single([]) );
    
    fprintf( '%i: calculate correlation\n', i )    
    c_it(i,:,:) = run_forward2_correlation_mex( mu, rho, G, spectrum, source_dist, rec, 0, [], single([]) );
    
    % [c_it(i,:,:), C] = run_forward2_correlation_mex( mu, rho, G, spectrum, source_dist, rec, 0, [], single([]) );    
    % parsave( ['../output/interferometry/C_' num2str(i) '.mat'], C )
    
    fprintf( '%i: done\n', i )
    
end
toc


% reorganize correlation vector
c_data = zeros(n_ref*n_rec,length(t));
for i = 1:n_ref
    c_data( (i-1)*n_rec + 1 : i*n_rec, :) = c_it(i,:,:);
end


% plot data
if( strcmp(make_plots,'yes') )
    figure
    plot_recordings(c_data,t,'vel','k-',true);
    legend('data')
end


% save array and data for inversion
save( sprintf('../output/interferometry/array_%i_ref.mat',n_ref), 'array', 'ref_stat')
save( sprintf('../output/interferometry/data_%i_ref_%i.mat',n_ref,n_basis_fct), 'c_data', 't')


% close matlabpool
if( ~strcmp(usr_par.cluster,'local') )
    delete(parobj)
end

