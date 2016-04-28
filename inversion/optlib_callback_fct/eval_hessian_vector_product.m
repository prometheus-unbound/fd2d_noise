% function [Hm, c_pert, c_ref] = eval_hessian_vector_product( m, dm, ModRandString, usr_par, mode )
function [Hm] = eval_hessian_vector_product( m, dm, ModRandString, usr_par, mode )
 

%- get configuration
[Lx, Lz, nx, nz, ~, ~, order, model_type, source_type, n_basis_fct] = input_parameters();
[~, ~, ~, ~, dx, dz] = define_computational_domain(Lx,Lz,nx,nz);

%- redirect optimization variable x and initialize kernel structures
m_parameters = map_m_to_parameters(m, usr_par);

% material parameters
mu = m_parameters(:,:,end);
[~,rho] = define_material_parameters( usr_par.config.nx, usr_par.config.nz, model_type );

% loading of spectra important for n_basis_fct=0, i.e. one map for each noise source
source_dist_ref = m_parameters(:,:,1:end-1);
[~, spectrum] = make_noise_source( source_type, n_basis_fct );

dm_parameters = map_m_to_parameters( m + dm, usr_par);
source_dist_pert = dm_parameters(:,:,1:end-1) - source_dist_ref;
dmu = dm_parameters(:,:,end) - mu;


%- loop over reference stations
t = usr_par.data.t;
nt = length(t);
n_ref = size( usr_par.network.ref_stat, 1 );
n_rec = size( usr_par.network.array, 1 ) - 1;
H_parameters = zeros(nx,nz,2);


parfor i = 1:n_ref
    
    
    %% each reference station will act as a source once
    src = usr_par.network.ref_stat(i,:);
    rec = usr_par.network.array( find( ~ismember(usr_par.network.array, src, 'rows') ) ,:);
    
    
    %% calculate Green function and perturbation
    [G_fft, G, ~, G_test] = run_forward1_green_mex( mu, rho, src, [], 5, [], single([]) );
    [dG_fft, dG, ~, dG_test] = run_forward1_green_mex( mu, rho, src, [], 5, dmu, G );

     
    %% calculate correlation and perturbation
    [c_ref, ~, C] = run_forward2_correlation_mex( mu, rho, G_fft, spectrum, source_dist_ref, rec, 1, [], single([]) );
    [c_pert1, ~, dC1] = run_forward2_correlation_mex( mu, rho, dG_fft, spectrum, source_dist_ref, rec, 1, dmu, C );
    [c_pert2, ~, dC2] = run_forward2_correlation_mex( mu, rho, G_fft, spectrum, source_dist_pert, rec, 1, [], single([]) );
    dC = dC1 + dC2;
    c_pert = c_pert1 + c_pert2;

    
    %% get corresponding data
    c_data = usr_par.data.c_data( (i-1)*n_rec + 1 : i*n_rec, : );
    
    
    %% calculate adjoint source time functions
    [~, adstf_1st_1_source] = make_adjoint_sources( c_ref, c_data, c_pert, t, usr_par.veldis, usr_par.measurement.source, src, rec, '1st' );
    [~, adstf_1st_1_structure] = make_adjoint_sources( c_ref, c_data, c_pert, t, usr_par.veldis, usr_par.measurement.structure, src, rec, '1st' );
    adstf_1st_1 = (1.0 - usr_par.kernel.weighting) * adstf_1st_1_source +  usr_par.kernel.weighting * adstf_1st_1_structure;
    
    [~, adstf_2nd_source] = make_adjoint_sources( c_ref, c_data, c_pert, t, usr_par.veldis, usr_par.measurement.source, src, rec, '2nd' );
    [~, adstf_2nd_structure] = make_adjoint_sources( c_ref, c_data, c_pert, t, usr_par.veldis, usr_par.measurement.structure, src, rec, '2nd' );
    adstf_2nd = (1.0 - usr_par.kernel.weighting) * adstf_2nd_source +  usr_par.kernel.weighting * adstf_2nd_structure;
    
    
    %% calculate adjoint wavefields and perturbations
    % u_dagger
    % dmu dmu - ds dmu
    [K_1, adstf_1st_2, u1] = run_noise_adjoint_mex( mu, rho, dC, complex(adstf_1st_1), rec, [], [], [], [], 13, [], single([]) );
    
    % p
    % dmu dmu
    [K_2_1, ~, p_1] = run_noise_adjoint_mex( mu, rho, dG, adstf_1st_2, rec, [], spectrum, source_dist_ref, [], 11, [], single([]) );
    
    % p
    % ds dmu
    [K_2_2, ~, p_2] = run_noise_adjoint_mex( mu, rho, G, adstf_1st_2, rec, [], spectrum, source_dist_pert, [], 11, [], single([]) );
    p = p_1 + p_2;
    
    % u_2_dagger
    % ds ds
    if( isempty( find( dmu, 1 ) ) )
        K_3_1 = run_noise_adjoint_mex( mu, rho, single([]), complex(adstf_2nd), rec, [], spectrum, [], G_fft, 10, [], single([]) );
    else
        K_3_1 = 0.0 * K_1;
    end

    % k
    % dmu dmu - ds dmu
    [K_3_2, adstf_pert] = run_noise_adjoint_mex( mu, rho, C, complex(adstf_2nd), rec, [], [], [], [], 12, dmu, u1 );
    
    % w und l
    % dmu dmu - ds dmu
    K_4 = run_noise_adjoint_mex( mu, rho, G, adstf_pert, rec, [], spectrum, source_dist_ref, [], 10, dmu, p );

    
    %% dmu ds
    % still a bit weird and not 100% correct - but don't know why
    K_5 = 0*K_4;
    K_6 = 0*K_4;
    
    if( isempty( find( source_dist_pert, 1 ) ) )
        n_noise_sources = size(spectrum,2);
        [~, n_sample, w_sample, dw, freq_samp] = input_interferometry();
        n_ftc = floor(nt/freq_samp);
        ifft_coeff = zeros(n_ftc,n_sample) + 1i*zeros(n_ftc,n_sample);
        
        i_ftc = 1;
        for n=1:nt
            if( mod(n,freq_samp) == 0 )
                for k = 1:n_sample
                    ifft_coeff(i_ftc,k) = 1/sqrt(2*pi) * exp( 1i*w_sample(k)*t(n) ) * dw;
                end
                i_ftc = i_ftc + 1;
            end
        end
        
        i_ftc = 1;
        n_zero = find(t==0);
        for n = 1:n_zero-1
            
            if( mod(n,freq_samp) == 0 && t(n) < 0.0 )
                
                T5 = zeros(nx,nz) + 1i*zeros(nx,nz);
                T6 = zeros(nx,nz) + 1i*zeros(nx,nz);
                for ns = 1:n_noise_sources
                    
                    for k = 1:n_sample
                        T5 = T5 + spectrum(k,ns) * conj(adstf_pert(:,:,k)) * ifft_coeff(i_ftc,k);
                        T6 = T6 + spectrum(k,ns) * conj(adstf_1st_2(:,:,k)) * ifft_coeff(i_ftc,k);
                    end
                    
                end
                
                i_ftc = i_ftc + 1;
                
                K_5(:,:,1) = K_5(:,:,1) + real(T5) .* G_test(:,:,end-n+1);
                K_6(:,:,1) = K_6(:,:,1) + real(T6) .* dG_test(:,:,end-n+1);
                
            end
        end
        
    end
    
    
    
    %% build up kernel
    H_parameters_i = K_1 + K_2_1 + K_2_2 + K_3_1 + K_3_2 + K_4 + K_5 + K_6;


    %% sum up kernels
    H_parameters = H_parameters + H_parameters_i;

    
end


%- map H_parameters to H_m
Hm = map_Hparameters_to_Hm( H_parameters, usr_par );


end
