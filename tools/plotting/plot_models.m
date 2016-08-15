% plot_models( structure, source, array, clim )
%
% structure
% source
% array (optional)
% clim: colormap limits (optional)
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
    width = width / 1000;
    array = array / 1000;
    
    
    % colormaps
    cm = cbrewer('div','RdBu',120,'PCHIP');
    cm_source = cm(50:120,:);

    
    % open figure, set size, etc.
    fig = figure;
    if( isempty( structure ) || isempty( source ) )
        set(fig,'units','normalized','position',[0.1 0.3 0.3 0.5])
    else
        set(fig,'units','normalized','position',[0.1 0.3 0.6 0.5])
        orient landscape
    end
    handle = [];
    legend_string = [];
    
    
    
    %% concerning structure
    if( ~isempty( structure ) )
        
        if( ~isempty( source ) )
            ax1 = subplot(1,2,1);
        else
            ax1 = gca;
            set(ax1,'position',[0.17 0.18 0.705 0.745]);
        end
        
        set(ax1,'FontSize',18);
        hold on
        
        % plot noise source
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
            handle(end+1,:) = plot3( array(:,1), array(:,2), level(1) + 0*array(:,2), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8 );
            legend_string{end+1} = 'array';
        end
        
        % include legend
        legend(handle,legend_string,'position',[0.44 0.06 0.1 0.05])
        
        % color stuff
        colormap(ax1,cm)
        if( clim(1)~=0 || clim(2)~=0 )
            caxis([clim(1) clim(2)]);
        end
        cb = colorbar;
        ylabel(cb,'[m/s]')
        m = max(max( abs(structure-4000) ));
        caxis([4000-m-double(m<1)*50 4000+m+double(m<1)*50])
        clabels = get(cb,'YTick');
        set(cb,'YTick',[clabels(1) clabels( ceil(length(clabels)/2) ) clabels(end)])
        
        % labels
        xlabels = [0 200 400];
        ylabels = [0 200 400];
        set(ax1, 'XTick', xlabels);
        set(ax1, 'YTick', ylabels);
        xlim([0 Lx])
        ylim([0 Lz])
        xlabel('x [km]')
        ylabel('z [km]')
        title('structure','FontSize',22)
        
        % layout
        shading interp
        grid on
        box on
        set(ax1,'LineWidth',2)
        axis square
        
        drawnow
        
    end
    
    
    
    %% concerning source
    if( ~isempty( source ) )
        
        if( ~isempty( structure ) )
            ax2 = subplot(1,2,2);
        else
            ax2 = gca;
            set(ax2,'position',[0.17 0.18 0.705 0.745]);
        end
        
        set(ax2,'FontSize',18);
        hold on
        
        % plot noise source
        mesh(X, Z, source' )
        
        % plot absorbing boundaries
        level = [1.1*max(max(source)), 1.1*max(max(source))];
        handle(end+1,:) = plot3([width,Lx-width],[width,width],level,'k--');
        plot3([width,Lx-width],[Lz-width,Lz-width],level,'k--')
        plot3([width,width],[width,Lz-width],level,'k--')
        plot3([Lx-width,Lx-width],[width,Lz-width],level,'k--')
        legend_string{end+1} = 'absorbing boundaries';
        
        % plot array if given
        if( ~isempty(array) )
            handle(end+1,:) = plot3( array(:,1), array(:,2), level(1) + 0*array(:,2), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8 );
            legend_string{end+1} = 'array';
        end
        
        % include legend
        if( isempty(structure) )
            legend(handle,legend_string,'position',[0.44 0.06 0.1 0.05])
        end
        
        % color stuff
        colormap(ax2,cm_source)
        if( clim(3)~=0 || clim(4)~=0 )
            caxis([clim(3) clim(4)]);
        end
        cb = colorbar;
        ylabel(cb,'[kg^2 m^{-2} s^{-2}]')        
        clabels = get(cb,'YTick');
        set(cb,'YTick',[clabels(1) clabels( ceil(length(clabels)/2) ) clabels(end)])
        
        % labels
        xlabels = [0 200 400];
        ylabels = [0 200 400];
        set(ax2, 'XTick', xlabels);
        xlim([0 Lx])
        ylim([0 Lz])
        xlabel('x [km]')
        if( isempty( structure ) )
            ylabel('z [km]')
            set(ax2, 'YTick', ylabels);
        else
            set(ax2, 'YTick', []);
        end
        
        title('source','FontSize',22)
        
        % layout
        shading interp
        grid on
        box on
        set(ax2,'LineWidth',2)
        axis square
        
        drawnow
        
    end    
    
    
end

