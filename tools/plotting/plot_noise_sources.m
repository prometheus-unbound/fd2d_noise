% plot_noise_sources(noise_source_distribution,array,cm_psd,clim,overlay)
%
% noise_source_distribution
% array (optional)
% cm_psd: colormap (optional)
% clim: colormap limits (optional)
% overlay: 'yes' or 'no', only of use for n_basis_fct == 0 (optional)
%
% optional: leave empty if not wanted

function [clim] = plot_noise_sources(noise_source_distribution,array,cm_psd,clim,overlay)

    % get configuration
    [f_sample,n_sample] = input_interferometry();
    [Lx,Lz,nx,nz,~,~,~,model_type,~,n_basis_fct] = input_parameters();
    [X,Z] = define_computational_domain(Lx,Lz,nx,nz);
    [width] = absorb_specs();
    X = X/1000; Lx = Lx/1000;
    Z = Z/1000; Lz = Lz/1000;
    width = width/1000;
    
    % if no colormap is given, load standard one
    if( isempty(cm_psd) )
        % load('cm_psd')
        cm = cbrewer('div','RdBu',120,'PCHIP');
        cm_psd = cm(52:120,:);
    end
 
    
    % open figure, set size, etc.
    fig1 = figure;
    % set(fig1,'units','normalized','position',[.1 -.3 0.5 1.1])
    set(gca,'FontSize',12);
    hold on
    
    
    % plot for n_basis_fct == 0, i.e. only one map
    if( n_basis_fct == 0 )
        
        if( nargin < 5 || isempty(overlay) )
            overlay = 'no';
        end
        
        if(strcmp(overlay,'yes'))

            [mu,~] = define_material_parameters(nx,nz,model_type);
            pcolor(X,Z,(mu-4.8e10)'/max(max(abs(mu-4.8e10))))

            dist = pcolor(X, Z, sum(noise_source_distribution,3)' / max(max( sum(noise_source_distribution,3) )) );
            alpha(dist,0.5)
            cb = colorbar;
            ylabel(cb,'normalized for overlay')
            caxis([-1.0 1.0])

        else
            
            if( size(noise_source_distribution,3) > 1 )
                pcolor(X, Z, sum(noise_source_distribution,3)'-1 )
            else
                pcolor(X, Z, sum(noise_source_distribution,3)' )
            end
            cb = colorbar;
            ylabel(cb,'relative source strength')

        end

        if( ~isempty(clim) )
            caxis(clim)
        else
            caxis([0.0 7.0])
            clim = get(gca,'CLim');
        end
        
%         plot([width,Lx-width],[width,width],'k--')
%         plot([width,Lx-width],[Lz-width,Lz-width],'k--')
%         plot([width,width],[width,Lz-width],'k--')
%         plot([Lx-width,Lx-width],[width,Lz-width],'k--')

        
    % plot maps for different frequency maps
    else

        fudge_factor = 10;

        int_limits = integration_limits(n_sample,n_basis_fct);
        for ib = 1:n_basis_fct

            h(ib) = mesh(X, Z, ib*fudge_factor + noise_source_distribution(:,:,ib)' );
            set(h(ib),'CData',noise_source_distribution(:,:,ib)');

            text(1e3,1e3, ib*fudge_factor + 1 + noise_source_distribution(1,1,ib) , sprintf('%5.3f - %5.3f Hz',f_sample(int_limits(ib,1)),f_sample(int_limits(ib,2))) )
            level = [ib*fudge_factor + 0.1 + noise_source_distribution(1,1,ib), ib*fudge_factor + 0.1 + noise_source_distribution(1,1,ib)];
            plot3([width,Lx-width],[width,width],level,'k--')
            plot3([width,Lx-width],[Lz-width,Lz-width],level,'k--')
            plot3([width,width],[width,Lz-width],level,'k--')
            plot3([Lx-width,Lx-width],[width,Lz-width],level,'k--')

        end

        colorbar
        
        if( ~isempty(clim) )
            caxis(clim)
        else
            clim = get(gca,'CLim');
        end

        view([22 4])
        % view([0 90])
        zlim([0 fudge_factor*(n_basis_fct+1)])
        set(gca,'ZTick',[])
        zlabel('frequency bands')

    end

    % plot array if given
    if( ~isempty(array) )
        plot(array(:,1),array(:,2),'ko')
    end
    
    xlabels = [0 1000 2000];
    ylabels = [0 1000 2000];
    set(gca, 'XTick', xlabels);
    set(gca, 'YTick', ylabels);
    
    shading interp   
    colormap(cm_psd);
    grid on
    box on
    axis square
    title('power-spectral density distribution');
    xlabel('x [km]')
    ylabel('z [km]')
    xlim([0 Lx])
    ylim([0 Lz])

    
end

