
%==========================================================================
% plot_kernel( gradient, usr_par )
%
% input:
%--------
% gradient
% usr_par
%
%==========================================================================


function plot_kernel(gradient, usr_par)


    %- configuration ------------------------------------------------------
    [Lx, Lz, nx, nz] = input_parameters();
    [X, Z] = define_computational_domain(Lx, Lz, nx, nz);
    [width, absorb_left, absorb_right, absorb_top, absorb_bottom] = absorb_specs();


    %- convert to km ------------------------------------------------------
    X = X / 1000; Z = Z / 1000;
    Lx = Lx / 1000; Lz = Lz / 1000;
    width = width / 1000;


    %- open figure, set size, etc. ----------------------------------------
    fig = figure;
    set(fig, 'units', 'normalized', 'position', [0.1, 0.3, 0.3, 0.5], 'Color', [1 1 1])
    ax1 = gca;
    set(ax1, 'FontSize', 12);
    hold on

    xlabel(ax1, 'x [km]')
    ylabel(ax1, 'z [km]')
    title(ax1, [usr_par.type, ' kernel'], 'FontSize', 14)

    handle = [];
    legend_string = [];


    %- plot gradient ------------------------------------------------------
    m = max(max(max(abs(gradient))));
    if (strcmp(usr_par.type, 'source'))
        pcolor(ax1, X, Z, gradient(:,:, 2)');
        caxis(ax1, [- 0.6 * m, 0.6 * m]);
    else
        pcolor(ax1, X, Z, gradient(:,:, 1)');
        caxis(ax1, [-0.1 * m, 0.1 * m]);
    end


    %- colormap and colorbar ----------------------------------------------
    cm = cbrewer('div', 'RdBu', 120);
    colormap(cm)
    cb = colorbar;
    clabels = get(cb, 'YTick');
    set(cb, 'YTick', [clabels(1), clabels(ceil(length(clabels) / 2)), clabels(end)])
    
    if (strcmp(usr_par.type, 'source'))
        ylabel(cb, '[kg^{-2} m^{2} s^{2}]')
    else
        ylabel(cb, '[N^{-1} m^{2}]')
    end


    %- plot absorbing boundaries ------------------------------------------
    % level = [100 * m, 100 * m];
    if(absorb_bottom); handle(1,:) = plot(ax1, [absorb_left*width, Lx - absorb_right*width], [width, width], 'k--'); end
        if(absorb_top); handle(1,:) = plot(ax1, [absorb_left*width, Lx - absorb_right*width], [Lz - width, Lz - width], 'k--'); end
        if(absorb_left); handle(1,:) = plot(ax1, [width, width], [absorb_bottom*width, Lz - absorb_top*width], 'k--'); end
        if(absorb_right); handle(1,:) = plot(ax1, [Lx - width, Lx - width], [absorb_bottom*width, Lz - absorb_top*width], 'k--'); end
    
    if( absorb_left || absorb_right || absorb_top || absorb_bottom );
        legend_string{end + 1} = 'absorbing boundaries';
    end


    %- plot array ---------------------------------------------------------
    handle(end + 1,:) = plot(ax1, usr_par.network.array(:, 1) / 1000, usr_par.network.array(:, 2) / 1000, ... % 1.1 * m + 0 * usr_par.network.array(:, 2), ...
        'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    
    
    %- legend -------------------------------------------------------------
    legend_string{end + 1} = 'array';
    legend(ax1, handle, legend_string)


    %- layout -------------------------------------------------------------
    xlabels = get(ax1, 'XTick');
    ylabels = get(ax1, 'YTick');
    set(ax1, 'XTick', [xlabels(1), xlabels(ceil(length(xlabels) / 2)), xlabels(end)])
    set(ax1, 'YTick', [ylabels(1), ylabels(ceil(length(ylabels) / 2)), ylabels(end)])
    shading(ax1, 'interp')
    xlim(ax1, [0, Lx])
    ylim(ax1, [0, Lz])
    grid(ax1, 'on')
    box(ax1, 'on')
    set(ax1, 'LineWidth', 2)
    axis(ax1, 'image')

    drawnow


end


