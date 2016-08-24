
clear all

%==========================================================================
% check adjoint source time functions
%==========================================================================

% load([fd2d_path() 'output' filesep 'array_nref-1.mat'])
% u = load([fd2d_path() 'output' filesep 'correlations_nref-1_model-1_source-homogeneous.mat']);
% u0 = load([fd2d_path() 'output' filesep 'correlations_nref-1_model-1_source-gaussian.mat']);
% 
% 
% %- choose du randomly -----------------------------------------------------
% du = 2 * (rand(1,length(u.t)) - 0.5);
% 
% %- choose one reference station and one receiver --------------------------
% i_ref = 1;
% rec_id = 1;
% 
% src = ref_stat(i_ref,:);
% rec = array(~ismember(array,src,'rows'), :);
% index = (i_ref - 1) * size(rec, 1) + 1 : i_ref * size(rec, 1);
% index = index(rec_id);
% rec = rec(rec_id,:);
% 
% 
% %- specify type of measurement and units of correlation -------------------
% usr_par.measurement.type = 'waveform_difference';
% % usr_par.measurement.type = 'log_amplitude_ratio';
% % usr_par.measurement.type = 'amplitude_difference';
% usr_par.measurement.mode = 'auto';
% 
% 
% %- run check --------------------------------------------------------------
% [dcheck, dcheck_struct] = optlib_check_adjoint_stf( u.correlations(index,:), u0.correlations(index,:), du, u.t, src, rec, -13, -2, 0.1, usr_par );


%==========================================================================
% wait till user starts gradient test 
%==========================================================================

% keyboard
% clear all


%==========================================================================
% check gradient
%==========================================================================

%- get basic configuration ------------------------------------------------
[Lx,Lz,nx,nz,dt,nt,order,model_type,source_type,n_basis_fct] = input_parameters();


%- set user parameters ----------------------------------------------------
usr_par.network = load([fd2d_path() 'output' filesep 'array_nref-1_testing.mat']);
usr_par.data = load([fd2d_path() 'output' filesep 'correlations_nref-1_testing.mat']);
usr_par.type = 'structure';
usr_par.measurement.type = 'waveform_difference';
usr_par.measurement.mode = 'auto';
usr_par.smoothing.sigma = [1e-3 1e-3];
usr_par.verbose = false;


%- set up initial model ---------------------------------------------------
noise_source = make_noise_source('no');
structure = define_material_parameters('no');

m_parameters = zeros(nx, nz, 2);
m_parameters(:,:,1) = structure.mu;
m_parameters(:,:,2) = noise_source.distribution;

m = map_parameters_to_m( m_parameters, usr_par );


%- include rho and spectrum into usr_par ----------------------------------
usr_par.structure.rho = structure.rho;
usr_par.noise_source.spectrum = noise_source.spectrum;


%- design test vector -----------------------------------------------------
dm = rand( numel(m), 1 );
dm( 1:nx*nz, 1 ) = 1e9 * dm( 1:nx*nz, 1 );


%- make sure that test vector is consistent with setup --------------------
if( strcmp( usr_par.type, 'structure' ) && ~isempty(find( dm(end-nx*nz+1:end, 1, 1), 1 )) )
    dm( nx*nz + 1 : end, 1 ) = 0 * dm( nx*nz + 1 : end, 1 );
elseif( strcmp( usr_par.type, 'source' ) && ~isempty(find( dm(1:end-nx*nz, 1), 1 )) )
    dm( 1:nx*nz, 1 ) = 0 * dm( 1:nx*nz, 1 );
end


%- set test vector in absorbing boundaries to zero ------------------------
[absbound] = reshape( init_absbound() , [], 1 );
for i = 1:size(dm)/(nx*nz)
   dm( ((i-1)*nx*nz+1):i*nx*nz, 1 ) = double( absbound == 1 ) .* dm( ((i-1)*nx*nz+1):i*nx*nz, 1 ); 
end


%- make sure that all necessary fields are set ----------------------------
[usr_par] = usr_par_init_default_parameters(usr_par);


%- run check --------------------------------------------------------------
[dcheck, dcheck_struct] = optlib_check_derivative( m, reshape(dm,[],1), -8, 0, 1, usr_par );

