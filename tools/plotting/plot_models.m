% plot_noise_sources(noise_source_distribution,array,cm_psd,clim,overlay)
%
% m_parameter
% array (optional)
% clim: colormap limits (optional)
% overlay: 'yes' or 'no', only of use for n_basis_fct == 0 (optional)
%
% optional: leave empty if not wanted

function plot_models( structure, source, array, clim )


    % configuration
    [Lx,Lz,nx,nz] = input_parameters();
    [X,Z] = define_computational_domain(Lx,Lz,nx,nz);
    [width] = absorb_specs();
    
    % convert to km
    X = X / 1000;  Lx = Lx / 1000;
    Z = Z / 1000;  Lz = Lz / 1000;
    width = width/1000;
    array = array / 1000;
    
    
    % colormaps
    cm = cbrewer('div','RdBu',120,'PCHIP');
    cm_source = cm(50:120,:);

    
    % open figure, set size, etc.
    fig1 = figure;
    set(fig1,'units','normalized','position',[.1 .3 0.8 0.6])
    orient landscape
    handle = [];
    legend_string = [];
    
    
    %% concerning source
    if( ~isempty( source ) )
        
        if( ~isempty( structure ) )
            ax1 = subplot(1,2,2);
        else
            ax1 = gca;
            set(ax1,'position',[0.35 0.11 0.3347 0.815]);
        end
        
        set(ax1,'FontSize',18);
        hold on
        
        mesh(X, Z, source' )
        
        % plot absorbing boundaries
        level = [1.1*max(max(source)), 1.1*max(max(source))];
        if( ~any(any( source<0 )) )
            handle(end+1,:) = plot3([width,Lx-width],[width,width],level,'k--');
            plot3([width,Lx-width],[Lz-width,Lz-width],level,'k--')
            plot3([width,width],[width,Lz-width],level,'k--')
            plot3([Lx-width,Lx-width],[width,Lz-width],level,'k--')
            legend_string{end+1} = 'absorbing boundaries';
        end
        
        % plot array if given
        if( ~isempty(array) )
            handle(end+1,:) = plot3( array(:,1), array(:,2), level(1) + 0*array(:,2), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8 );
            legend_string{end+1} = 'array';
        end
        
        legend(handle,legend_string,'position',[0.44 0.06 0.1 0.05])
        cb = colorbar;
        ylabel(cb,'[kg^2 m^{-2} s^{-2}]')
        
        if( any(any( source<0 )) )
            colormap(ax1,cm)
        else
            colormap(ax1,cm_source)
        end
        
        if( clim(3)~=0 || clim(4)~=0 )
            caxis([clim(3) clim(4)]);
        elseif( any(any( source<0 )) )
            m = max(max(abs(source)));
            caxis([-0.6*m 0.6*m])
        end
        
        clabels = get(cb,'YTick');
        set(cb,'YTick',[clabels(1) clabels( ceil(length(clabels)/2) ) clabels(end)])
        
        xlabels = [0 200 400];
        ylabels = [0 200 400];
        set(ax1, 'XTick', xlabels);
        set(ax1, 'YTick', ylabels);
        xlim([0 Lx])
        ylim([0 Lz])
        xlabel('x [km]')
        ylabel('z [km]')
        title('source','FontSize',24)
        
        shading interp
        grid on
        box on
        set(ax1,'LineWidth',2)
        axis square
        
    end
    
    
    if( isempty( structure ) )
        return
    end
    
    
    
    %% concerning structure
    if( ~isempty( source ) )
        ax2 = subplot(1,2,1);
    else
        ax2 = gca;
        set(ax2,'position',[0.35 0.11 0.3347 0.815]);
    end
    
    set(ax2,'FontSize',18);
    hold on
    
    mesh(X, Z, structure' )
    
    % plot absorbing boundaries
    level = [1.1*max(max(structure)), 1.1*max(max(structure))];
    handle(end+1,:) = plot3([width,Lx-width],[width,width],level,'k--');
    plot3([width,Lx-width],[Lz-width,Lz-width],level,'k--')
    plot3([width,width],[width,Lz-width],level,'k--')
    plot3([Lx-width,Lx-width],[width,Lz-width],level,'k--')
    legend_string{end+1} = 'absorbing boundaries';
    
    % plot array if given
    if( ~isempty(array) )
        handle(end+1,:)  = plot3( array(:,1), array(:,2), level(1) + 0*array(:,2), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8 );
        legend_string{end+1} = 'array';
    end
    
    if( isempty( source ) )
        legend(handle,legend_string,'position',[0.44 0.06 0.1 0.05])
    end
    
    cb = colorbar;
    ylabel(cb,'[m/s]')
    colormap(ax2,cm);
    
    if( clim(1)~=0 || clim(2)~=0 )
        caxis([clim(1) clim(2)]);
        set(cb,'YTick',[clim(1) clim(2)])
    else
        
    end
    
    clabels = get(cb,'YTick');
    set(cb,'YTick',[clabels(1) clabels( ceil(length(clabels)/2) ) clabels(end)])
    
    xlabels = [0 200 400];
    ylabels = [0 200 400];
    set(ax2, 'XTick', xlabels);
    if( ~isempty( source ) )
        set(ax2, 'YTick', []);
    else
        set(ax2, 'YTick', ylabels);
    end
    xlim([0 Lx])
    ylim([0 Lz])
    xlabel('x [km]')
    title('structure','FontSize',24)
    
    shading interp   
    grid on
    box on
    set(ax2,'LineWidth',2)
    axis square
    
    drawnow
    
    
end

