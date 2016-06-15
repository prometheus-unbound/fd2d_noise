
function [Hm] = eval_hessian_vector_product_source( m, dm, ModRandString, usr_par )
 

if( ~strcmp( usr_par.type, 'source' ) )
    error('\nonly use that function for usr_par.type = "source"\n')
end


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
    dm_generic( 1:end-nx*nz, 1 ) = dm;
else
    dm_generic = dm;
end

dm_parameters = map_m_to_parameters( m + dm_generic, usr_par);
source_dist_pert = dm_parameters(:,:,1:end-1) - source_dist_ref;


%- loop over reference stations
t = usr_par.data.t;
n_ref = size( usr_par.network.ref_stat, 1 );
n_rec = size( usr_par.network.array, 1 ) - 1;
H_parameters = 0 * m_parameters;


parfor i = 1:n_ref
    
    
    %% each reference station will act as a source once
    src = usr_par.network.ref_stat(i,:);
    rec = usr_par.network.array( find( ~ismember(usr_par.network.array, src, 'rows') ) ,:);
    
    
    %% get corresponding data
    c_data = usr_par.data.c_data( (i-1)*n_rec + 1 : i*n_rec, : );
    
    
    %% calculate Green function and perturbation
    if( exist(['../output/interferometry/G_fft_' num2str(i) '.mat'], 'file') )
        G_fft = parload( ['../output/interferometry/G_fft_' num2str(i) '.mat'] );
    else
        [G_fft] = run_forward1_green_mex(mu, rho, src, [], 0, [], single([]) );
        parsave( ['../output/interferometry/G_fft_' num2str(i) '.mat'], G_fft )
        if( strcmp( usr_par.verbose, 'yes') ) fprintf( '%i: Green function -- done\n', i ); end
    end
    
     
    %% calculate correlation and perturbation
    if( ~strcmp( usr_par.measurement.source, 'waveform_difference' ) || ~strcmp( usr_par.measurement.structure, 'waveform_difference' ) )
        
        if( exist(['../output/interferometry/c_ref_' num2str(i) '.mat'], 'file') )
            c_ref = parload( ['../output/interferometry/c_ref_' num2str(i) '.mat'] );
        else
            c_ref = run_forward2_correlation_mex( mu, rho, G_fft, spectrum, source_dist_ref, rec, 0, [], single([]) );
            parsave( ['../output/interferometry/c_ref_' num2str(i) '.mat'], c_ref )
            if( strcmp( usr_par.verbose, 'yes') ) fprintf( '%i: correlation -- done\n', i ); end
        end
        
    else
        c_ref = 0.0 * c_data;
    end
    
    
    % perturbation of correlation following from source perturbation
    c_pert = run_forward2_correlation_mex( mu, rho, G_fft, spectrum, source_dist_pert, rec, 0, [], single([]) );
    if( strcmp( usr_par.verbose, 'yes') ) fprintf( '%i: correlation perturbation -- done\n', i ); end
    
    
    %% calculate adjoint source time functions
    [~, adstf_2nd_source] = make_adjoint_sources( c_ref, c_data, c_pert, t, usr_par.veldis, usr_par.measurement.source, src, rec, '2nd' );
    [~, adstf_2nd_structure] = make_adjoint_sources( c_ref, c_data, c_pert, t, usr_par.veldis, usr_par.measurement.structure, src, rec, '2nd' );
    adstf_2nd = (1.0 - usr_par.kernel.weighting) * adstf_2nd_source +  usr_par.kernel.weighting * adstf_2nd_structure;
    
    
    %% calculate adjoint wavefields and perturbations    
    % BLACK -- k
    H_parameters_i = run_noise_adjoint_mex( mu, rho, single([]), complex(adstf_2nd), rec, [], spectrum, [], G_fft, 10, [], single([]) );
    if( strcmp( usr_par.verbose, 'yes') ) fprintf( '%i: BLACK -- done\n', i ); end
    

    %% sum up kernels
    H_parameters = H_parameters + H_parameters_i;

    
end


%- map H_parameters to H_m
Hm = map_Hparameters_to_Hm( H_parameters, usr_par );


end
