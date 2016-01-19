
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% usr_par.type = 'source';
usr_par.type = 'structure';

usr_par.measurement = 'waveform_difference';
% 'log_amplitude_ratio';
% 'amplitude_difference';
% 'waveform_difference';
% 'cc_time_shift';

usr_par.network = load('../output/interferometry/array_1_ref.mat');
usr_par.data = load('../output/interferometry/data_1_ref_0.mat');

usr_par.data_independent = 'yes';

usr_par.use_mex = 'no';

[usr_par] = usr_par_init_default_parameters_lbfgs([]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up model
[~,~,nx,nz,dt,nt,~,model_type,~,n_basis_fct] = input_parameters();

% get source and material
[source_dist, spectrum] = make_noise_source();
[mu, rho] = define_material_parameters(nx,nz,model_type);

% specify output behaviour
output_specs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- loop over reference stations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = -(nt-1)*dt:dt:(nt-1)*dt;
nt = length(t);

n_ref = size( usr_par.network.ref_stat, 1 );
n_rec = size( usr_par.network.array, 1 ) - 1;

c_it = zeros( n_ref, n_rec, nt );

j = 0;

if(n_basis_fct == 0)
    n_basis_fct = 1;
end
grad_parameters = zeros( nx, nz, n_basis_fct);


for i = 1:n_ref
    
    
    % each reference station will act as a source once
    src = usr_par.network.ref_stat(i,:);
    rec = usr_par.network.array( find( ~ismember(usr_par.network.array, src, 'rows') ) ,:);
    
    
    % load or calculate Green function
    if( strcmp( usr_par.type, 'source') && exist(['../output/interferometry/G_2_' num2str(i) '.mat'], 'file') )
        
        G_2 = load_G_2( ['../output/interferometry/G_2_' num2str(i) '.mat'] );
        
    else
        
        if( strcmp( usr_par.type, 'source') )
            
            if( strcmp( usr_par.use_mex, 'yes') )
                [G_2] = run_forward1_green_mex(mu, rho, src, 0);
            else
                [G_2] = run_forward1_green(mu, rho, src, 0);
            end            
            parsave( ['../output/interferometry/G_2_' num2str(i) '.mat'], G_2 )

        else
            
            if( strcmp( usr_par.use_mex, 'yes') )
                [G_2, G_2_dxu_time, G_2_dzu_time] = run_forward1_green_mex(mu, rho, src, 1);
            else
                [G_2, G_2_dxu_time, G_2_dzu_time] = run_forward1_green(mu, rho, src, 1);
            end
            
        end
    end
    
    
    % calculate correlation
    if( strcmp( usr_par.type, 'source') )
        
        if( strcmp( usr_par.use_mex, 'yes') )
            [c_it(i,:,:)] = run_forward2_correlation_mex( mu, rho, G_2, spectrum, source_dist, rec, 0, usr_par.debug.df );
        else
            [c_it(i,:,:)] = run_forward2_correlation( mu, rho, G_2, spectrum, source_dist, rec, 0, usr_par.debug.df );
        end
        
    elseif( strcmp( usr_par.type, 'structure') )
        
        if( strcmp( usr_par.use_mex, 'yes') )
            [c_it(i,:,:), ~, C_2_dxu_time, C_2_dzu_time] = run_forward2_correlation_mex( mu, rho, G_2, spectrum, source_dist, rec, 1, 0 );
        else
            [c_it(i,:,:), ~, C_2_dxu_time, C_2_dzu_time] = run_forward2_correlation( mu, rho, G_2, spectrum, source_dist, rec, 1, 0 );
        end
        
    end
    
    
    % filter correlations if wanted
    if( strcmp( usr_par.filter.apply_filter, 'yes') )
        c_data_iref = filter_correlations( usr_par.data.c_data( (i-1)*n_rec + 1 : i*n_rec, : ), t, usr_par.filter.f_min, usr_par.filter.f_max, 4 );
        c_it(i,:,:) = filter_correlations( reshape(c_it(i,:,:),[],nt), t, usr_par.filter.f_min, usr_par.filter.f_max, 4 );
               
    else
            
%         if( strcmp(usr_par.data_independent, 'yes' ))
            c_data_iref = 0 * c_it;
%         else
%             c_data_iref = usr_par.data.c_data( (i-1)*n_rec + 1 : i*n_rec, : );
%         end
        
    end
    
    
    % calculate misfit and adjoint source function
    [j_n, adstf] = make_adjoint_sources_inversion( reshape(c_it(i,:,:),[],nt), c_data_iref, t, usr_par.veldis, usr_par.measurement, src, rec );
    
    
    % calculate source kernel
    if( strcmp( usr_par.type, 'source') )
        
        if( strcmp( usr_par.use_mex, 'yes') )
            grad_parameters_i = run_noise_source_kernel_mex( mu, rho, G_2, spectrum, adstf, rec );
        else
            grad_parameters_i = run_noise_source_kernel( mu, rho, G_2, spectrum, adstf, rec );
        end
        
    % calculate structure kernel
    elseif( strcmp( usr_par.type, 'structure') )
        
        if( strcmp( usr_par.use_mex, 'yes') )
            
            % first run
            [grad_mu_i_1, adjoint_state_1] = run_noise_structure_kernel_mex( mu, rho, C_2_dxu_time, C_2_dzu_time, complex(adstf), rec, spectrum, source_dist );
            
            % second run
            [grad_mu_i_2] = run_noise_structure_kernel_mex( mu, rho, G_2_dxu_time, G_2_dzu_time, adjoint_state_1, rec, spectrum, source_dist );
            
        else
            
            % first run
            [grad_mu_i_1, adjoint_state_1] = run_noise_structure_kernel( mu, rho, C_2_dxu_time, C_2_dzu_time, complex(adstf), rec, spectrum, source_dist );
            
            % second run
            [grad_mu_i_2] = run_noise_structure_kernel( mu, rho, G_2_dxu_time, G_2_dzu_time, adjoint_state_1, rec, spectrum, source_dist );
            
        end
        
        % sum up both contributions
        grad_parameters_i = grad_mu_i_1 + grad_mu_i_2;
        
    end
    
    
    % sum up kernels
    grad_parameters = grad_parameters + grad_parameters_i;
    
    
    % sum up misfits
    j = j + sum(j_n);
    
    
end


%- reorganize correlation vector
c_all = zeros( n_ref*n_rec, length(t) );
for i = 1:n_ref
    c_all( (i-1)*n_rec + 1 : i*n_rec, :) = c_it(i,:,:);
end


figure
mesh(grad_parameters')
cm = cbrewer('div','RdBu',120,'PCHIP');
colormap(cm)
c = max( max( grad_parameters ));
caxis( [-0.7*c 0.7*c] )
colorbar

% % plot data and synthetics
% if (strcmp(make_plots,'yes'))
%     figure
%     hold on
%     if( strcmp(data_independent,'yes') )
%         h(1,:) = plot_recordings(c_all,t,'vel','r-',true);
%         legend(h,'uniform')
%     else
%         h(1,:) = plot_recordings(c_data,t,'vel','k-',true);
%         h(2,:) = plot_recordings(c_uniform,t,'vel','r-',true);
%         legend(h,'data','uniform')
%     end
%     drawnow
% end


