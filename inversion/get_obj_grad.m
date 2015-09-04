
function [f, g, c_all] = get_obj_grad(x)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    measurement = 1;
    % 1 = 'log_amplitude_ratio';
    % 2 = 'amplitude_difference';
    % 3 = 'waveform_difference';
    % 4 = 'cc_time_shift';
    
    % do measurement on displacement or velocity correlations 
    % (for NOW: use 'dis')
    veldis = 'dis';
    
    % filter correlations for inversion; if yes, specify f_min and f_max
    apply_filter = 'yes';
    f_min = 1/7 - 0.005;
    f_max = 1/7 + 0.005;
    
    % load array with reference stations and data
    load('../output/interferometry/array_16_ref.mat');
    load('../output/interferometry/data_16_ref_0_uniform_2gaussian_homogeneous.mat');
    
    % specify percentile for clipping of kernel      ( = 0 to turn it off )
    percentile = 0;                                 
    
    % desing gaussian filter for smoothing of kernel ( = 0 to turn it off )
    % myfilter = 0;
    myfilter = fspecial('gaussian',[75 75], 30);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate misfit and gradient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    global type;
    global v0;

    [~,~,nx,nz,dt,nt,~,model_type] = input_parameters();
    f_sample = input_interferometry(); 

    %- redirect optimization variable x and initialize kernel structures   
    if( strcmp(type,'source') )
        
        % loading of spectra important for n_basis_fct=0, i.e. one map for each noise source
        source_dist = x;
        [~,spectrum] = make_noise_source();
        
        
        % load structure that is assumed for the source inversion
        [mu,~] = define_material_parameters(nx,nz,model_type); 
        
        % allocate kernel
        K_all = zeros( nx, nz, length(f_sample) );
        
        
    elseif( strcmp(type,'structure') )
        
        % get source that is assumed for the structure inversion
        [source_dist,spectrum] = make_noise_source();
        
        
        % need density for reparameterization from mu to v kernel
        [~,rho] = define_material_parameters( nx, nz, model_type);
        
        % get mu from v, which is our optimization variable
        mu = v0^2 * reshape( rho, [], 1) .* (1+x).^2;
        
        % allocate kernel
        K_all = zeros(nx, nz);
        
    end
    
    
    
    %- loop over reference stations
    f = 0;
    t = -(nt-1)*dt:dt:(nt-1)*dt;
    nt = length(t);
    n_ref = size(ref_stat,1);
    n_rec = size(array,1)-1;
    c_it = zeros(n_ref,n_rec,nt);
    
    
    parfor i = 1:n_ref
        
        
        % each reference station will act as a source once
        src = ref_stat(i,:);
        rec = array( find( ~ismember(array, src, 'rows') ) ,:);
        
        
        % load or calculate Green function
        if( strcmp(type,'source') && exist(['../output/interferometry/G_2_' num2str(i) '.mat'], 'file') )
            G_2 = load_G_2( ['../output/interferometry/G_2_' num2str(i) '.mat'] );
            
        else
            [G_2] = run_forward_green_fast_mex(mu, src);
            
            if( strcmp(type,'source') )
                parsave( ['../output/interferometry/G_2_' num2str(i) '.mat'], G_2 )
            end
        end
        
        
        % calculate correlation
        if( strcmp(type,'source') )
            [c_it(i,:,:), ~] = run_forward_correlation_fast_mex( G_2, source_dist, spectrum, mu, rec, 0 );
        elseif( strcmp(type,'structure') )
            [c_it(i,:,:), ~, C_2_dxv, C_2_dzv] = run_forward_correlation_fast_mex( G_2, source_dist, spectrum, mu, rec, 1 );
        end
        
        
        % filter correlations if wanted
        if( strcmp(apply_filter,'yes') )
            c_data_iref = filter_correlations( c_data( (i-1)*n_rec + 1 : i*n_rec, : ), t, f_min, f_max );
            c_it(i,:,:) = filter_correlations( reshape(c_it(i,:,:),[],nt), t, f_min, f_max );
        end
        
        
        % calculate misfit and adjoint source function
        switch measurement
            case 1
                [f_n,adstf] = make_adjoint_sources_inversion( reshape(c_it(i,:,:),[],nt), c_data_iref, t, veldis, 'log_amplitude_ratio', src, rec );
            case 2
                [f_n,adstf] = make_adjoint_sources_inversion( reshape(c_it(i,:,:),[],nt), c_data_iref, t, veldis, 'amplitude_difference', src, rec );
            case 3
                [f_n,adstf] = make_adjoint_sources_inversion( reshape(c_it(i,:,:),[],nt), c_data_iref, t, veldis, 'waveform_difference', src, rec );
            case 4
                [f_n,adstf] = make_adjoint_sources_inversion( reshape(c_it(i,:,:),[],nt), c_data_iref, t, veldis, 'cc_time_shift', src, rec );
            otherwise
                error('\nspecify correct measurement!\n\n')
        end
        
        
        % calculate source kernel
        if( strcmp(type,'source') )                
            [~,~,K_i] = run_noise_source_kernel_fast_mex( G_2, mu, adstf, rec );
        
        % calculate structure kernel    
        elseif( strcmp(type,'structure') )                
            [~,~,K_i] = run_noise_mu_kernel_fast_mex( C_2_dxv, C_2_dzv, mu, adstf, rec );
        
        end
        
        
        % sum up kernels
        K_all = K_all + K_i;
        
        % sum up misfits
        f = f + f_n;
        
        
    end
    
    
    fprintf('misfit: %f\n',f)
    
    
    %- reorganize correlation vector
    c_all = zeros( n_ref*n_rec, length(t) );
    for i = 1:n_ref
        c_all( (i-1)*n_rec + 1 : i*n_rec, :) = c_it(i,:,:);
    end
    
    
    %- sum source kernel (all frequencies or in frequency bands)
    if( strcmp(type,'source') )
        K_all = sum_source_kernel(K_all);
    end
    
    
    %- kernel treatment, i.e. clipping/smoothing
    K_all = treat_kernel( K_all, type, myfilter, percentile );
    
    
    %- reshape kernel to column vector 
    if( strcmp(type,'source') )
        g = reshape( K_all, [], 1 );
    elseif( strcmp(type,'structure') )
        g = 2 * v0^2 * reshape( rho .* K_all, [], 1 ) .* (1+x);
    end

      
end


