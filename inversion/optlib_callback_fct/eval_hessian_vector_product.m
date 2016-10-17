
function [Hm] = eval_hessian_vector_product( m, dm, ModRandString, usr_par )
 

%- get configuration
[~, ~, nx, nz, ~, ~, ~, model_type, source_type, n_basis_fct] = input_parameters();

%- redirect optimization variable x and initialize kernel structures
m_parameters = map_m_to_parameters(m, usr_par);

% material parameters
mu = m_parameters(:,:,end);
[~,rho] = define_material_parameters( usr_par.config.nx, usr_par.config.nz, model_type );

% loading of spectra important for n_basis_fct=0, i.e. one map for each noise source
source_dist_ref = m_parameters(:,:,1:end-1);
[~, spectrum] = make_noise_source( source_type, n_basis_fct );

if( numel(m_parameters) ~= numel(dm) )

    dm_generic = 0 * m;
    
    if( strcmp( usr_par.type, 'source' ) )
        dm_generic( 1:end-nx*nz, 1 ) = dm;
    elseif( strcmp( usr_par.type, 'structure' ) )
        dm_generic( end-nx*nz+1:end, 1 ) = dm;
    end
    
else
    dm_generic = dm;
end

dm_parameters = map_m_to_parameters( m + dm_generic, usr_par);
source_dist_pert = dm_parameters(:,:,1:end-1) - source_dist_ref;
dmu = dm_parameters(:,:,end) - mu;


%- loop over reference stations
t = usr_par.data.t;
n_ref = size( usr_par.network.ref_stat, 1 );
n_rec = size( usr_par.network.array, 1 ) - 1;
H_parameters = 0 * m_parameters;


for i = 1:n_ref
    
    
    %% each reference station will act as a source once
    src = usr_par.network.ref_stat(i,:);
    rec = usr_par.network.array( find( ~ismember(usr_par.network.array, src, 'rows') ) ,:);
    
    
    %% get corresponding data
    c_data = usr_par.data.c_data( (i-1)*n_rec + 1 : i*n_rec, : );
    
    
    %% calculate Green function and perturbation
    [G_fft, G] = run_forward1_green_mex( mu, rho, src, [], 1, [], single([]) );
    
    % perturbation of Green function following from dmu
    if( ~isempty( find( dmu, 1 ) ) )
        [dG_fft, dG] = run_forward1_green_mex( mu, rho, src, [], 5, dmu, G );
        if( strcmp( usr_par.verbose, 'yes') ) fprintf( '%i: dGreen dmu -- done\n', i ); end
    else
        dG_fft = 0.0;
    end
    
    if( strcmp( usr_par.verbose, 'yes') ) fprintf( '%i: Green function -- done\n', i ); end

     
    %% calculate correlation and perturbation
    [c_ref, C] = run_forward2_correlation_mex( mu, rho, G_fft, spectrum, source_dist_ref, rec, 1, [], single([]) );
    
    % perturbation of correlation following from dmu
    if( ~isempty( find( dmu, 1 ) ) )
        [c_pert1, dC] = run_forward2_correlation_mex( mu, rho, dG_fft, spectrum, source_dist_ref, rec, 1, dmu, single(C) );
        if( strcmp( usr_par.verbose, 'yes') ) fprintf( '%i: dcorrelation dmu -- done\n', i ); end
    else
        c_pert1 = 0.0;
        dC = 0.0;
    end
    
    % perturbation of correlation following from source perturbation
    if( ~isempty( find( source_dist_pert, 1 ) ) )
        [c_pert2, dC2] = run_forward2_correlation_mex( mu, rho, G_fft, spectrum, source_dist_pert, rec, 1, [], single([]) );
        if( strcmp( usr_par.verbose, 'yes') ) fprintf( '%i: dcorrelation dsource -- done\n', i ); end
    else
        c_pert2 = 0.0;
        dC2 = 0.0;
    end    
    
    dC = dC + dC2;
    dC2 = [];
    c_pert = c_pert1 + c_pert2;
    
    if( strcmp( usr_par.verbose, 'yes') ) fprintf( '%i: correlation -- done\n', i ); end
    
    
    %% calculate adjoint source time functions
    [~, adstf_1st_1_source] = make_adjoint_sources( c_ref, c_data, c_pert, t, usr_par.veldis, usr_par.measurement.source, src, rec, '1st' );
    [~, adstf_1st_1_structure] = make_adjoint_sources( c_ref, c_data, c_pert, t, usr_par.veldis, usr_par.measurement.structure, src, rec, '1st' );
    adstf_1st_1 = (1.0 - usr_par.kernel.weighting) * adstf_1st_1_source +  usr_par.kernel.weighting * adstf_1st_1_structure;
    
    [~, adstf_2nd_source] = make_adjoint_sources( c_ref, c_data, c_pert, t, usr_par.veldis, usr_par.measurement.source, src, rec, '2nd' );
    [~, adstf_2nd_structure] = make_adjoint_sources( c_ref, c_data, c_pert, t, usr_par.veldis, usr_par.measurement.structure, src, rec, '2nd' );
    adstf_2nd = (1.0 - usr_par.kernel.weighting) * adstf_2nd_source +  usr_par.kernel.weighting * adstf_2nd_structure;
    
    
    %% calculate adjoint wavefields and perturbations
    K_1 = 0.0;  K_2 = 0.0;  K_4 = 0.0;
    
    %% YELLOW -- u_dagger
    % dmu dmu - ds dmu - dmu ds
    if( ~strcmp( usr_par.type, 'source' ) )
        [K_1, uT, u_dag] = run_noise_adjoint_mex( mu, rho, single(dC), complex(adstf_1st_1), rec, [], spectrum, [], dG_fft, 13, [], single([]) );
    else
        u_dag = 0.0;
    end
    
    if( strcmp( usr_par.type, 'structure' ) )
        K_1(:,:,1:end-1) = 0 * K_1(:,:,1:end-1);
    end
    
    dC = [];  dG_fft = [];
    if( strcmp( usr_par.verbose, 'yes') ) fprintf( '%i: YELLOW -- done\n', i ); end
    
    
    %% GRAY -- p
    % dmu dmu
    if( ~strcmp( usr_par.type, 'source' ) && ~isempty( find( dmu, 1 ) ) )
        [K_2, ~, p] = run_noise_adjoint_mex( mu, rho, single(dG), complex(uT), rec, [], spectrum, source_dist_ref, complex([]), 11, [], single([]) );
    else
        p = 0.0;
    end
    
    dG = [];
    if( strcmp( usr_par.verbose, 'yes') ) fprintf( '%i: GRAY 1 -- done\n', i ); end
    
    
    %% GRAY -- p
    % ds dmu
    if( strcmp( usr_par.type, 'joint' ) && ~isempty( find( source_dist_pert, 1 ) ) )
        [K_2_add, ~, p_add] = run_noise_adjoint_mex( mu, rho, single(G), complex(uT), rec, [], spectrum, source_dist_pert, complex([]), 11, [], single([]) );
        K_2 = K_2 + K_2_add;
        p = p + p_add;
        p_add = [];
    end
    
    uT = [];
    if( strcmp( usr_par.verbose, 'yes') ) fprintf( '%i: GRAY 2 -- done\n', i ); end
    
    
    %% BLACK -- k
    % dmu dmu - ds dmu - dmu ds - ds ds
    [K_3, kT] = run_noise_adjoint_mex( mu, rho, single(C), complex(adstf_2nd), rec, [], spectrum, [], G_fft, 12, dmu, single(u_dag) );
    
    if( strcmp( usr_par.type, 'source' ) )
        K_3(:,:,end) = 0 * K_3(:,:,end);
    elseif( strcmp( usr_par.type, 'structure') )
        K_3(:,:,1:end-1) = 0 * K_3(:,:,1:end-1);
    end
    
    u_dag = [];  C = [];  G_fft = [];
    if( strcmp( usr_par.verbose, 'yes') ) fprintf( '%i: BLACK -- done\n', i ); end
    
    
    %% RED -- w und l
    % dmu dmu - ds dmu
    if( ~strcmp( usr_par.type, 'source') )
        K_4 = run_noise_adjoint_mex( mu, rho, single(G), complex(kT), rec, [], spectrum, source_dist_ref, complex([]), 10, dmu, single(p) );
    end
    
    if( strcmp( usr_par.verbose, 'yes') ) fprintf( '%i: RED -- done\n', i ); end
    
    
    %% build up kernel
    H_parameters_i = K_1 + K_2 + K_3 + K_4;


    %% sum up kernels
    H_parameters = H_parameters + H_parameters_i;

    
end


%- map H_parameters to H_m
Hm = map_Hparameters_to_Hm( H_parameters, usr_par );


%- add term from regularization - Skpye-Moji saying: WHAAAAATTT ?????
% [absbound] = reshape( init_absbound() , [], 1 );
% Hm( 1:end - nx*nz, 1 ) = Hm( 1:end - nx*nz, 1 ) + double( absbound == 1 ) .* usr_par.regularization.alpha * usr_par.regularization.weighting;
% Hm( nx*nz+1:end, 1 ) = Hm( nx*nz+1:end, 1 ) + double( absbound == 1 ) .* usr_par.regularization.beta * usr_par.regularization.weighting;


end
