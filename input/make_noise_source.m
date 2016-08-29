
function [ noise_source_distribution, noise_spectrum ] = make_noise_source( source_type, n_basis_fct, make_plots )

    % define make_plots if not specified
    if( nargin < 3 )
        make_plots = 'no';
    end

    
    % get configuration
    [f_sample, n_sample] = input_interferometry();
    [Lx, Lz, nx, nz, ~, ~, ~, model_type] = input_parameters();
    [X, Z] = define_computational_domain(Lx,Lz,nx,nz);  
    
    if( Lx == 2.0e6 || Lz == 2.0e6 )
        size = 'big';
    else
        size = 'small';
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % user input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    n_noise_sources = 1;

    
    %- specify spectrum ---------------------------------------------------
    
    % two overlapping spectra
    f_peak = [1/11, 1/7];
    bandwidth = [0.035, 0.025];
    strength = [0.3, 1];

    % narrow non-overlapping spectra
    % f_peak = [1/15, 1/7];
    % bandwidth = [0.01, 0.01];
    % strength = [1, 1];

    
    
    %- different source types --------------------------------------------- 
    
    % location and width of a gaussian 'blob'
    if(strcmp(source_type,'gaussian'))
        
        % small setup, 2 sources 
        if( strcmp(size,'small') )
            
            x_sourcem = [1.0e5];
            z_sourcem = [1.6e5];
            sourcearea_width = [4e4];
            magnitude = [6.0];
            
        % large setup, 2 sources left of the array
        elseif( strcmp(size,'big') )
            
            % standard
            x_sourcem = [0.6e6, 0.5e6, 1.3e6];
            z_sourcem = [0.8e6, 1.5e6, 1.6e6];
            sourcearea_width = [2.0e5, 1.5e5, 1.8e5];
            magnitude = [6.0, 5.0, 3.0];
            
            % x_sourcem = [0.4e6];
            % z_sourcem = [1.0e6];
            % sourcearea_width = [0.3e5];
            % magnitude = [6.0];
            
            % laura
            % x_sourcem = 1.3e6;
            % z_sourcem = 0.8e6;
            % sourcearea_width = 1.5e5;
            % magnitude = 6.0;
                      
            % ring test
            % x_sourcem = 0.5e6;
            % z_sourcem = 0.8e6;
            % sourcearea_width = 1.5e5;
            % magnitude = 6.0;

        end


    % ring of sources
    elseif(strcmp(source_type,'ring'))
        
        x_source_r = 1.0e6;
        z_source_r = 1.0e6;
        radius = 6.8e5;
        thickness = 1e5;
        angle_cover = 60.0;
        taper_width = 20.0;
        taper_strength = 100;
    
        
    % picture translated sources
    elseif(strcmp(source_type,'picture'))
        
        filename = 'source.png';
        
    end
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % source spectrum
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    noise_spectrum = zeros(length(f_sample),n_noise_sources);
    
    for ns = 1:n_noise_sources
        noise_spectrum(:,ns) = strength(ns) * exp( -(abs(f_sample)-f_peak(ns)).^2 / bandwidth(ns)^2 );       
    end
    
%     if( strcmp(source_type, 'homogeneous') )
%         noise_spectrum = sum(noise_spectrum,2);
%         n_noise_sources = 1;
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % geographic power-spectral density distribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % homogeneous noise source distribution
    noise_source_distribution = ones(nx,nz,n_noise_sources);
    
    
    % gaussian blob
    if(strcmp(source_type,'gaussian'))
        
        for ns = 1:n_noise_sources
            noise_source_distribution(:,:,ns) = noise_source_distribution(:,:,ns) ...
                + magnitude(ns) * ( exp( -( (X-x_sourcem(ns)).^2 + (Z-z_sourcem(ns)).^2 ) / (sourcearea_width(ns))^2 ) )';
        end
        
%         for ns = 1:3
%             noise_source_distribution(:,:,1) = noise_source_distribution(:,:,1) ...
%                 + magnitude(ns) * ( exp( -( (X-x_sourcem(ns)).^2 + (Z-z_sourcem(ns)).^2 ) / (sourcearea_width(ns))^2 ) )';
%         end

%         for ns = 1
%             noise_source_distribution(:,:,1) = noise_source_distribution(:,:,1) ...
%                 + magnitude(ns) * ( exp( -(X-0.5e6).^2  / 1.8e5^2 ) .* exp( -(Z-1.0e6).^2 / 5.0e5^2 ) )';
%         end

        
    % ring of sources with taper
    elseif( strcmp(source_type,'ring') )
        
        R = ( (X-x_source_r).^2 + (Z-z_source_r).^2 ).^(1/2);
        angle = atan( abs( X-x_source_r ) ./ abs( Z-z_source_r ) ) *180/pi;
        
        [k,l] = find(isnan(angle));
        angle(k,l) = 0;
        
        if( angle_cover == 90 )
            noise_source_distribution(:,:,1) = (exp( -abs( R-radius ).^2/9e8 ) .* double(R > (radius-thickness/2) & R < (radius+thickness/2) ) );
        else
            noise_source_distribution(:,:,1) = (exp( -abs( R-radius ).^2/9e8 ) .* double(R > (radius-thickness/2) & R < (radius+thickness/2) ) );
            noise_source_distribution(:,:,1) = noise_source_distribution(:,:,1) ...
                + exp( -abs( R-radius ).^2/9e8 ) .* ( 5*exp( -(angle-(angle_cover-taper_width)).^2/(taper_strength) ...
                .* double( angle>angle_cover-taper_width & angle<angle_cover ) ) ...
                .* double( R > (radius-thickness/2) & R < (radius+thickness/2) & angle <= angle_cover ) );
        end
            
            
    % translate picture to source
    elseif(strcmp(source_type,'picture'))
        
        A = imread(filename);
        noise_source_distribution(:,:,1) = flipud( abs((double(A(:,:))-255) / max(max(abs(double(A)-255)))) )';
        
        
    end
    
    
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % translate spectrum and distribution to source appropriate for
    % run_forward
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if( n_basis_fct ~= 0 )
        [noise_source_distribution] = general_source( noise_spectrum, noise_source_distribution, n_basis_fct );
        noise_spectrum = ones(n_sample,1);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot noise source configuration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ( strcmp(make_plots,'yes') )
        
        if ( n_basis_fct == 0 )
            
            cmap = hsv(n_noise_sources);
            for ns = 1:n_noise_sources
                
                if(ns==1)
                    figure
                    set(gca,'FontSize',12)
                    hold on
                    grid on
                    cstring = [];
                end
                
                plot(f_sample,noise_spectrum(:,ns),'Color',cmap(ns,:))
                cstring{end+1} = ['source ' num2str(ns)];
                xlabel('frequency [Hz]');
                legend(cstring)
                
                xlim([f_sample(1) f_sample(end)])
                
            end
            
        end
        
        
        % load ../output/interferometry/array_28_ref.mat
        % load ../output/interferometry/array_1_ref_hessian_big.mat
        % load ../output/interferometry/array_16_ref.mat
        % load ../output/interferometry/array_16_ref_small.mat
        % load ../output/interferometry/array_16_ref_center2.mat
        array = [];
        
        if( n_basis_fct == 0 )
            m_parameters = zeros( nx, nz, 2);
        else
            m_parameters = zeros( nx, nz, n_basis_fct+1 );
        end
        
        m_parameters(:,:,1:end-1) = noise_source_distribution;
        m_parameters(:,:,end) = define_material_parameters( nx, nz, model_type );
        
        usr_par.network = [];
        usr_par.data = [];
        usr_par.config.n_basis_fct = n_basis_fct;
        
        % usr_par.kernel.imfilter.source = fspecial('gaussian',[75 75], 30);
        % usr_par.kernel.imfilter.source = fspecial('gaussian',[40 40], 20);
        usr_par.kernel.sigma.source = [1 1];
        usr_par.kernel.sigma.structure = usr_par.kernel.sigma.source;
        
        [usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);
        m_parameters = map_m_to_parameters( map_parameters_to_m(m_parameters, usr_par ) , usr_par );
        
        % cm = cbrewer('div','RdBu',120,'PCHIP');
        % plot_models( m_parameters, n_basis_fct, array, [-0.2 1.1 3.8e3 4.2e3], 'no', 'no', cm );       
        plot_models( m_parameters, n_basis_fct, array, [0 7 3.8e3 4.2e3], 'no', 'no' );       
        
    end
    
    
end

