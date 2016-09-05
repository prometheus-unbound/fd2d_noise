
%==========================================================================
% plot_models( structure, source, array, clim )
%
% input:
%--------
% structure
% source
% array (optional)
% clim: colormap limits (optional)
%
%==========================================================================

function plot_models(structure, source, array, clim)


    %- configuration ------------------------------------------------------
    [Lx, Lz, nx, nz] = input_parameters();
    [X, Z] = define_computational_domain(Lx, Lz, nx, nz);
    [width, absorb_left, absorb_right, absorb_top, absorb_bottom] = absorb_specs();


    %- convert to km ------------------------------------------------------
    X = X / 1000; Lx = Lx / 1000;
    Z = Z / 1000; Lz = Lz / 1000;
    width = width / 1000;
    array = array / 1000;


    %- colormaps ----------------------------------------------------------
    cm = cbrewer('div', 'RdBu', 120);
    cm_source = cm(50:120,:);


    %- open figure, set size, etc. ----------------------------------------
    fig = figure;
    if (isempty(structure) || isempty(source))
        set(fig, 'units', 'normalized', 'position', [0.1, 0.3, 0.3, 0.5])
    else
        set(fig, 'units', 'normalized', 'position', [0.1, 0.3, 0.6, 0.5])
        orient landscape
    end


    %======================================================================
    % concerning structure
    %======================================================================

    if (~isempty(structure))

        handle = [];
        legend_string = [];
        
        if (~isempty(source))
            ax1 = subplot(1, 2, 1);
        else
            ax1 = gca;
            set(ax1, 'position', [0.17, 0.18, 0.705, 0.745]);
        end

        set(ax1, 'FontSize', 12);
        hold(ax1, 'on')


        %- plot noise source ----------------------------------------------
        mesh(ax1, X, Z, structure')


        %- plot absorbing boundaries --------------------------------------
        level = [1.1 * max(max(structure)), 1.1 * max(max(structure))];
        if(absorb_bottom); handle(1,:) = plot3(ax1, [absorb_left*width, Lx - absorb_right*width], [width, width], level, 'k--'); end
        if(absorb_top); handle(1,:) = plot3(ax1, [absorb_left*width, Lx - absorb_right*width], [Lz - width, Lz - width], level, 'k--'); end
        if(absorb_left); handle(1,:) = plot3(ax1, [width, width], [absorb_bottom*width, Lz - absorb_top*width], level, 'k--'); end
        if(absorb_right); handle(1,:) = plot3(ax1, [Lx - width, Lx - width], [absorb_bottom*width, Lz - absorb_top*width], level, 'k--'); end
        
        if( absorb_left || absorb_right || absorb_top || absorb_bottom );
            legend_string{end + 1} = 'absorbing boundaries';
        end


        %- plot array if given --------------------------------------------
        if (~isempty(array))
            handle(end + 1,:) = plot3(ax1, array(:, 1), array(:, 2), level(1) + 0 * array(:, 2), ...
                'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
            legend_string{end + 1} = 'array';
        end


        %- include legend -------------------------------------------------
        legend(ax1, handle, legend_string, 'position', [0.44, 0.06, 0.1, 0.05])


        %- color stuff ----------------------------------------------------
        colormap(ax1, cm)
        if (clim(1) ~= 0 || clim(2) ~= 0)
            caxis(ax1, [clim(1), clim(2)]);
        end
        cb = colorbar;
        ylabel(cb, '[m/s]')
        m = max(max(abs(structure - 4000)));
        caxis(ax1, [4000 - m - double(m < 1) * 50, 4000 + m + double(m < 1) * 50])
        clabels = get(cb, 'YTick');
        set(cb, 'YTick', [clabels(1), clabels(ceil(length(clabels) / 2)), clabels(end)])


        %- labels ---------------------------------------------------------
        xlabels = get(ax1, 'XTick');
        ylabels = get(ax1, 'YTick');
        set(ax1, 'XTick', [xlabels(1), xlabels(ceil(length(xlabels) / 2)), xlabels(end)])
        set(ax1, 'YTick', [ylabels(1), ylabels(ceil(length(ylabels) / 2)), ylabels(end)])
        xlim(ax1, [0, Lx])
        ylim(ax1, [0, Lz])
        xlabel(ax1, 'x [km]')
        ylabel(ax1, 'z [km]')
        title(ax1, 'structure', 'FontSize', 14)


        %- layout ---------------------------------------------------------
        shading(ax1, 'interp')
        grid(ax1, 'on')
        box(ax1, 'on')
        set(ax1, 'LineWidth', 2)
        axis(ax1, 'image')


    end


    %======================================================================
    % concerning source
    %======================================================================
    
    if (~isempty(source))
        
        handle = [];
        legend_string = [];
        
        if (~isempty(structure))
            ax2 = subplot(1, 2, 2);
        else
            ax2 = gca;
            set(ax2, 'position', [0.17, 0.18, 0.705, 0.745]);
        end

        set(ax2, 'FontSize', 12);
        hold(ax2, 'on')


        %- plot noise source ----------------------------------------------
        mesh(ax2, X, Z, source')


        %- plot absorbing boundaries --------------------------------------
        level = [1.1 * max(max(source)), 1.1 * max(max(source))];
        if(absorb_bottom); handle(1,:) = plot3(ax2, [absorb_left*width, Lx - absorb_right*width], [width, width], level, 'k--'); end
        if(absorb_top); handle(1,:) = plot3(ax2, [absorb_left*width, Lx - absorb_right*width], [Lz - width, Lz - width], level, 'k--'); end
        if(absorb_left); handle(1,:) = plot3(ax2, [width, width], [absorb_bottom*width, Lz - absorb_top*width], level, 'k--'); end
        if(absorb_right); handle(1,:) = plot3(ax2, [Lx - width, Lx - width], [absorb_bottom*width, Lz - absorb_top*width], level, 'k--'); end
        
        if( absorb_left || absorb_right || absorb_top || absorb_bottom );
            legend_string{end + 1} = 'absorbing boundaries';
        end


        %- plot array if given --------------------------------------------
        if (~isempty(array))
            handle(end + 1,:) = plot3(ax2, array(:, 1), array(:, 2), level(1) + 0 * array(:, 2), ...
                'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
            legend_string{end + 1} = 'array';
        end


        %- include legend -------------------------------------------------
        if (isempty(structure))
            legend(ax2, handle, legend_string, 'position', [0.44, 0.06, 0.1, 0.05])
        end


        %- color stuff ----------------------------------------------------
        colormap(ax2, cm_source)
        if (clim(3) ~= 0 || clim(4) ~= 0)
            caxis(ax2, [clim(3), clim(4)]);
        end
        cb = colorbar;
        ylabel(cb, '[kg^2 m^{-2} s^{-2}]')
        clabels = get(cb, 'YTick');
        set(cb, 'YTick', [clabels(1), clabels(ceil(length(clabels) / 2)), clabels(end)])


        %- labels ---------------------------------------------------------
        xlabels = get(ax2, 'XTick');
        set(ax2, 'XTick', [xlabels(1), xlabels(ceil(length(xlabels) / 2)), xlabels(end)])
                        
        xlim(ax2, [0, Lx])
        ylim(ax2, [0, Lz])
        xlabel(ax2, 'x [km]')
        if (isempty(structure))
            ylabel(ax2, 'z [km]')
            ylabels = get(ax2, 'YTick');
            set(ax2, 'YTick', [ylabels(1), ylabels(ceil(length(ylabels) / 2)), ylabels(end)])
        else
            set(ax2, 'YTick', []);
        end

        title(ax2, 'source', 'FontSize', 14)


        %- layout ---------------------------------------------------------
        shading(ax2, 'interp')
        grid(ax2, 'on')
        box(ax2, 'on')
        set(ax2, 'LineWidth', 2)
        axis(ax2, 'image')


    end
    
    
    drawnow


end


