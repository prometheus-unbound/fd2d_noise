
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
    % Lx = Lx / 1000; Lz = Lz / 1000;
    % width = width / 1000;


    %- open figure, set size, etc. ----------------------------------------
    fig = figure;
    set(fig, 'units', 'normalized', 'position', [0.1, 0.3, 0.3, 0.5])
    set(gca, 'FontSize', 16);
    hold on

    xlabel('x [km]')
    ylabel('z [km]')
    set(gca, 'XTick', [0, 200, 400])
    set(gca, 'YTick', [0, 200, 400])
    title([usr_par.type, ' kernel'], 'FontSize', 20)

    handle = [];
    legend_string = [];


    %- plot gradient ------------------------------------------------------
    m = max(max(max(abs(gradient))));
    if (strcmp(usr_par.type, 'source'))
        mesh(X, Z, gradient(:,:, 2)');
        caxis([- 0.6 * m, 0.6 * m]);
    else
        mesh(X, Z, gradient(:,:, 1)');
        caxis([- 0.1 * m, 0.1 * m]);
    end


    %- colormap and colorbar ----------------------------------------------
    cm = cbrewer('div', 'RdBu', 120, 'PCHIP');
    colormap(cm)
    cb = colorbar;
    clabels = get(cb, 'YTick');
    set(cb, 'YTick', [clabels(1), clabels(ceil(length(clabels) / 2)), clabels(end)])


    %- plot absorbing boundaries ------------------------------------------
    % level = [1.1 * m, 1.1 * m];
    % handle(end + 1,:) = plot3([width, Lx - width], [width, width], level, 'k--');
    % plot3([width, Lx - width], [Lz - width, Lz - width], level, 'k--')
    % plot3([width, width], [width, Lz - width], level, 'k--')
    % plot3([Lx - width, Lx - width], [width, Lz - width], level, 'k--')
    % legend_string{end + 1} = 'absorbing boundaries';


    %- legend -------------------------------------------------------------
    handle(end + 1,:) = plot3(usr_par.network.array(:, 1) / 1000, usr_par.network.array(:, 2) / 1000, 1.1 * m + 0 * usr_par.network.array(:, 2), ...
        'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    legend_string{end + 1} = 'array';
    legend(handle, legend_string)


    %- layout -------------------------------------------------------------
    shading interp
    grid on
    box on
    set(gca, 'LineWidth', 2)
    axis square

    drawnow


end


