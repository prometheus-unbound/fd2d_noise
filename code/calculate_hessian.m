
addpath(genpath('../'))
[~,~,nx,nz,~,~,~,model_type,source_type,n_basis_fct] = input_parameters();
usr_par.config.nx = nx; 
usr_par.config.nz = nz;
usr_par.config.n_basis_fct = n_basis_fct; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

usr_par.cluster = 'local';
% 'local';
% 'monch';
% 'euler';
% 'brutus';


usr_par.type = 'joint';
% 'source'
% 'structure'
% 'joint'


usr_par.use_mex = 'no';
% 'yes'
% 'no'


% define reference models for perturbations and mu_0 for structure
usr_par.initial.ref_source = 0;
usr_par.initial.ref_structure = 1;
usr_par.initial.mu_0 = 4.8e10;


% if wanted, load initial perturbations for source and/or structure
% initial_model = load('initial_models/source_random_0.10_loga_equal_homogeneous_regu_1em4_smooth_5e4.mat');
% usr_par.initial.source_dist = initial_model.model.m( 1 : end - nx*nz, 1 );
% initial_model = load('initial_models/structure_random_0.10_cc_equal_homogeneous_regu_1em2_smooth_5e4.mat');
% usr_par.initial.structure = initial_model.model.m( end - nx*nz + 1 : end, 1 );


usr_par.measurement.source = 'log_amplitude_ratio';
usr_par.measurement.structure = 'waveform_difference';
% 'log_amplitude_ratio';
% 'amplitude_difference';
% 'waveform_difference';
% 'cc_time_shift';


% define weighting for kernels
usr_par.kernel.weighting = 0.5;


% design gaussian filter for smoothing of kernel (set second input variable to [1,1] to turn it off)
% usr_par.kernel.imfilter.source = fspecial('gaussian',[75 75], 30);
% usr_par.kernel.imfilter.source = fspecial('gaussian',[1 1], 1);
% usr_par.kernel.imfilter.source = fspecial('gaussian',[40 40], 20);
% usr_par.kernel.imfilter.source = fspecial('gaussian',[1 1], 1);
% usr_par.kernel.imfilter.structure = usr_par.kernel.imfilter.source;

usr_par.kernel.sigma.source = [1e-3 1e-3];
usr_par.kernel.sigma.structure = usr_par.kernel.sigma.source;


% load array with reference stations and data
% usr_par.network = load('../output/interferometry/array_16_ref.mat');
% usr_par.data = load('../output/interferometry/data_16_ref_0_gaussian_random_0.10_0.8e10_nosmooth.mat');

usr_par.network = load('../output/interferometry/array_2_ref_testing.mat');
usr_par.data = load('../output/interferometry/data_2_ref_0_testing.mat');


% do measurement on displacement or velocity correlations (for NOW: use 'dis')
usr_par.veldis = 'dis';


% regularization j_total = dj/dm + alpha * ||m_source-m0_source||^2 + beta * ||m_structure-m0_structure||^2 
% (i.e. set alpha or beta to zero to turn it off)
usr_par.regularization.alpha = 0;
usr_par.regularization.beta = 0;
usr_par.regularization.weighting = weighting( nx, nz );


usr_par.verbose = 'yes';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% running the inversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( strcmp( usr_par.type, 'source') )
    
    if( usr_par.kernel.weighting ~= 0.0 )
        
        usr_par.kernel.weighting = 0.0;
        usr_par.regularization.beta = 0.0;
        
        fprintf('\nset usr_par.kernel.weighting to 0.0')
        fprintf('\nset usr_par.regularization.beta to 0.0\n')
        
    end
    
elseif( strcmp( usr_par.type, 'structure') )
    
    if( usr_par.kernel.weighting ~= 1.0 )
        
        usr_par.kernel.weighting = 1.0;
        usr_par.regularization.alpha = 0.0;
    
        fprintf('\nset usr_par.kernel.weighting to 1.0')
        fprintf('\nset usr_par.regularization.alpha to 0.0\n')
    end
    
end


% set necessary fields that might not have been set
[usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);


% set up initial model
if( usr_par.config.n_basis_fct == 0 )
    m_parameters = zeros( nx, nz, 2 );
else
    m_parameters = zeros( nx, nz, usr_par.config.n_basis_fct + 1 );
end

m_parameters(:,:,1:end-1) = make_noise_source( source_type, usr_par.config.n_basis_fct );
m_parameters(:,:,end) = define_material_parameters( nx, nz, model_type );


% convert to optimization variable
m0 = map_parameters_to_m( m_parameters, usr_par );


if( isfield( usr_par, 'initial') )
    if( isfield( usr_par.initial, 'source_dist') )
        m0( 1 : end - nx*nz, 1 ) = usr_par.initial.source_dist;
    end
    
    if( isfield( usr_par.initial, 'structure') )
        m0( end - nx*nz + 1 : end, 1 ) = usr_par.initial.structure;
    end
end


% initial model for regularization
usr_par.m0 = m0;


% set up test vector
dm = 0 * m_parameters;
% dm( 382:386, 293:297, 2 ) = 0.1;
dm( 25:30, 25:30, 2 ) = 0.1;


if( strcmp( usr_par.cluster, 'local' ) )
    
    figure
    [Lx,Lz] = input_parameters();
    [X,Z] = define_computational_domain(Lx,Lz,nx,nz);
    [width] = absorb_specs();
    offset = max( max( max( abs(dm) ) ) );
    
    hold on
    mesh( X, Z, dm(:,:,1)' )
    plot3( usr_par.network.array(:,1), usr_par.network.array(:,2), 0*usr_par.network.array(:,1) + offset, 'x' )
    plot3([width,Lx-width],[width,width],[offset,offset],'k--')
    plot3([width,Lx-width],[Lz-width,Lz-width],[offset,offset],'k--')
    plot3([width,width],[width,Lz-width],[offset,offset],'k--')
    plot3([Lx-width,Lx-width],[width,Lz-width],[offset,offset],'k--')
    axis square
    
%     return
    
end


% set test vectors in absorbing boundary region to zero
dm = reshape( dm, [], 1 );
[absbound] = reshape( init_absbound() , [], 1 );
for i = 1:size(dm)/(nx*nz)
   dm( ((i-1)*nx*nz+1):i*nx*nz, 1 ) = double( absbound == 1 ) .* dm( ((i-1)*nx*nz+1):i*nx*nz, 1 ); 
end


% tmp = reshape(dm,600,600,2);
% mesh( X, Z, tmp(:,:,2)' )
% return


% start matlabpool
if( ~strcmp( usr_par.cluster, 'local') )
    parobj = start_cluster( usr_par.cluster, '', size(usr_par.network.ref_stat,1));
end


% calculate Hessian vector product
tic
Hdm = eval_hessian_vector_product( m0, dm, optlib_generate_random_string(8), usr_par );
toc


% plot Hessian vector product
if( strcmp( usr_par.cluster, 'local' ) )
    
    fig666 = figure;%(666);
    set(fig666,'units','normalized','position',[.1 .5 0.5 0.4])
    
    [Lx,Lz] = input_parameters();
    [X,Z] = define_computational_domain(Lx,Lz,nx,nz);
    
    uiui = reshape( Hdm, nx, nz, 2 );
    angle = [0 90];
    
    subplot(1,2,1)
    mesh( X, Z, uiui(:,:,1)' )
    cm = cbrewer('div','RdBu',120,'PCHIP');
    colormap(cm);
    caxis([ -max(max(abs( uiui(:,:,1) )))  max(max(abs( uiui(:,:,1) )))]);
%     axis image
    colorbar
    view(angle)
    
    subplot(1,2,2)
    mesh( X, Z, uiui(:,:,end)' )
    cm = cbrewer('div','RdBu',120,'PCHIP');
    colormap(cm);
    caxis([ -max(max(abs( uiui(:,:,end) )))  max(max(abs( uiui(:,:,end) )))]);
%     axis image
    colorbar
    view(angle)
    
    drawnow
    
end


% save solution
save( '../output/hessian.mat', 'Hdm', 'usr_par', 'dm', '-v7.3' )


% close matlabpool
if( ~strcmp( usr_par.cluster, 'local' ) )
    delete(parobj)
end

