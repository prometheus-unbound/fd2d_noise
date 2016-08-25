
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
    [width] = absorb_specs();


    %- convert to km ------------------------------------------------------
    X = X / 1000; Z = Z / 1000;
    Lx = Lx / 1000; Lz = Lz / 1000;
    width = width / 1000;


    %- open figure, set size, etc. ----------------------------------------
    fig = figure;
    set(fig, 'units', 'normalized', 'position', [0.1, 0.3, 0.3, 0.5])
    ax1 = gca;
    set(ax1, 'FontSize', 12);
    hold on

    xlabel(ax1, 'x [km]')
    ylabel(ax1, 'z [km]')
    set(ax1, 'XTick', [0, 200, 400])
    set(ax1, 'YTick', [0, 200, 400])
    title(ax1, [usr_par.type, ' kernel'], 'FontSize', 14)

    handle = [];
    legend_string = [];


    %- plot gradient ------------------------------------------------------
    m = max(max(max(abs(gradient))));
    if (strcmp(usr_par.type, 'source'))
        mesh(ax1, X, Z, gradient(:,:, 2)');
        caxis(ax1, [- 0.6 * m, 0.6 * m]);
    else
        mesh(ax1, X, Z, gradient(:,:, 1)');
        caxis(ax1, [- 0.1 * m, 0.1 * m]);
    end


    %- colormap and colorbar ----------------------------------------------
    cm = cbrewer('div', 'RdBu', 120);
    colormap(cm)
    cb = colorbar;
    clabels = get(cb, 'YTick');
    set(cb, 'YTick', [clabels(1), clabels(ceil(length(clabels) / 2)), clabels(end)])


    %- plot absorbing boundaries ------------------------------------------
    level = [1.1 * m, 1.1 * m];
    handle(end + 1,:) = plot3(ax1, [width, Lx - width], [width, width], level, 'k--');
    plot3(ax1, [width, Lx - width], [Lz - width, Lz - width], level, 'k--')
    plot3(ax1, [width, width], [width, Lz - width], level, 'k--')
    plot3(ax1, [Lx - width, Lx - width], [width, Lz - width], level, 'k--')
    legend_string{end + 1} = 'absorbing boundaries';


    %- plot array ---------------------------------------------------------
    handle(end + 1,:) = plot3(ax1, usr_par.network.array(:, 1) / 1000, usr_par.network.array(:, 2) / 1000, 1.1 * m + 0 * usr_par.network.array(:, 2), ...
        'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    
    
    %- legend -------------------------------------------------------------
    legend_string{end + 1} = 'array';
    legend(ax1, handle, legend_string)


    %- layout -------------------------------------------------------------
    shading(ax1, 'interp')
    grid(ax1, 'on')
    box(ax1, 'on')
    set(ax1, 'LineWidth', 2)
    axis(ax1, 'square')

    drawnow


end


