
function [noise_source_distribution,noise_spectrum] = make_noise_source(make_plots)

    % define make_plots if not specified
    if( nargin < 1 )
        make_plots = 'no';
    end

    % get configuration
    [f_sample,n_sample] = input_interferometry();
    [Lx,Lz,nx,nz,~,~,~,model_type,source_type,n_basis_fct] = input_parameters();
    [X,Z] = define_computational_domain(Lx,Lz,nx,nz);
    [width] = absorb_specs();
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % user input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    n_noise_sources = 2;
    
    
    %- specify spectrum ---------------------------------------------------
    
    % two overlapping spectra
    f_peak = [1/11, 1/7];
    bandwidth = [0.035, 0.025];
    strength = [1, 0.7];

%     % narrow non-overlapping spectra
%     f_peak = [1/15, 1/7];
%     bandwidth = [0.01, 0.01];
%     strength = [1, 0.7];

%     % one simple source
%     f_peak = 0.125;
%     bandwidth = 0.03;
%     strength = 1;
    
    
    %- different source types --------------------------------------------- 
    
    % location and width of a gaussian 'blob'
    if(strcmp(source_type,'gaussian'))
        
%         % small setup, 2 sources 
%         x_sourcem = [2.0e5, 1.0e5];
%         z_sourcem = [2.0e5, 1.0e5];
%         sourcearea_width = [0.4e5, 0.4e5];
%         magnitude = [3.0, 2.0];
        
        % large setup, 2 sources left of the array
        x_sourcem = [0.5e6, 0.6e6];
        z_sourcem = [0.8e6, 1.3e6];
        sourcearea_width = [2.0e5, 1.5e5];
        magnitude = [6.0, 5.0];
        
%         % large steup, 1 source in the center of the array
%         x_sourcem = [1.25e6];
%         z_sourcem = [1.0e6];
%         sourcearea_width = [1.5e5];
%         magnitude = [3.0];

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
    
    cmap = hsv(6);
    for ns = 1:n_noise_sources
        noise_spectrum(:,ns) = strength(ns) * exp( -(abs(f_sample)-f_peak(ns)).^2 / bandwidth(ns)^2 );
        
        if ( strcmp(make_plots,'yes') )
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
            
            % figure
            % mesh(noise_source_distribution(:,:,ns))
        end
        
        
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
        [noise_source_distribution] = general_source( noise_spectrum, noise_source_distribution );
        noise_spectrum = ones(n_sample,1);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot noise distribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ( strcmp(make_plots,'yes') )
        
        if( n_basis_fct == 0 )
            
            overlay = 'no';
            
            figure
            set(gca,'FontSize',12);
            hold on
            
            if(strcmp(overlay,'yes'))
                
                [mu,~] = define_material_parameters(nx,nz,model_type);
                pcolor(X,Z,(mu-4.8e10)'/max(max(abs(mu-4.8e10))))
                
                dist = pcolor(X, Z, sum(noise_source_distribution,3)' / max(max( sum(noise_source_distribution,3) )) );
                alpha(dist,0.5)
                cm = cbrewer('div','RdBu',100,'PCHIP');
                colormap(cm);
                cb = colorbar;
                ylabel(cb,'normalized for overlay')
                caxis([-1.0 1.0])
                
            else

                pcolor(X, Z, sum(noise_source_distribution,3)' )
                load cm_psd
                colormap(cm_psd)
                colorbar
                
            end
            
            plot([width,Lx-width],[width,width],'k--')
            plot([width,Lx-width],[Lz-width,Lz-width],'k--')
            plot([width,width],[width,Lz-width],'k--')
            plot([Lx-width,Lx-width],[width,Lz-width],'k--')
                    
        else
            
            fudge_factor = 10;
            
            fig1 = figure;
            set(fig1,'units','normalized','position',[.1 -.3 0.5 1.1])
            set(gca,'FontSize',12);
            hold on
            
            int_limits = integration_limits(n_sample,n_basis_fct);
            
            for ib = 1:n_basis_fct
                
                noise_dist_basis = noise_source_distribution(:,:,ib);
                h(ib) = mesh(X, Z, ib*fudge_factor + noise_dist_basis' );
                set(h(ib),'CData',noise_dist_basis');
                
                text(1e3,1e3, ib*fudge_factor + 1 + noise_dist_basis(1,1) , sprintf('%5.2f - %5.2f Hz',f_sample(int_limits(ib,1)),f_sample(int_limits(ib,2))) )
                level = [ib*fudge_factor + 0.1 + noise_dist_basis(1,1), ib*fudge_factor + 0.1 + noise_dist_basis(1,1)];
                plot3([width,Lx-width],[width,width],level,'k--')
                plot3([width,Lx-width],[Lz-width,Lz-width],level,'k--')
                plot3([width,width],[width,Lz-width],level,'k--')
                plot3([Lx-width,Lx-width],[width,Lz-width],level,'k--')
                
            end
            
            load cm_psd
            colormap(cm_psd)
            colorbar
            
            view([22 4])
            zlim([0 fudge_factor*(n_basis_fct+1)])
            set(gca,'ZTick',[])
            zlabel('frequency bands')
            
        end
        
        % load('~/Desktop/runs/inversion/data/array_16_ref.mat')
        % plot(array(:,1),array(:,2),'ko')        
        
        shading interp
        grid on
        box on
        axis square
        title('power-spectral density distribution');
        xlabel('x [m]')
        ylabel('z [m]')
        xlim([0 Lx])
        ylim([0 Lz])
        
    end
    
end

