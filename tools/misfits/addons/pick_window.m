
%==========================================================================
% [left,right] = pick_window( u, u_0, t, measurement )
%==========================================================================


function [left, right] = pick_window(u, u_0, t, measurement)


    %- convert to velocity ------------------------------------------------
    % v = zeros( 1, length(t) );
    % v_0 = zeros( 1, length(t) );

    % v(1,1:end-1) = diff( u ) / (t(2)-t(1));
    % v_0(1,1:end-1) = diff( u_0 ) / (t(2)-t(1));


    %- filter displacement recordings -------------------------------------
    u = filter_seismogram(u, t, 0.02, 0.2, 4);
    u_0 = filter_seismogram(u_0, t, 0.02, 0.2, 4);


    %- open figure, set size, etc. ----------------------------------------
    fig = figure;
    set(fig, 'units', 'normalized', 'position', [0.1, 0.3, 0.6, 0.5])
    set(gca, 'FontSize', 14);
    hold on
    box on
    legend_string = [];


    %- plot seismograms ---------------------------------------------------
    plot(t, u, 'r')
    legend_string{end + 1} = 'synthetics';
    if( sum(u_0) ~= 0 )
        plot(t, u_0, 'k')
        legend_string{end + 1} = 'data';
    end


    % legend, labels, title, limits ---------------------------------------
    legend(legend_string)
    xlabel('time [s]')
    title('pick measurement window', 'FontSize', 16)
    xlim([t(1), t(end)])
    ylim([- 1.2 * max(max(abs(u)), max(abs(u_0))), 1.2 * max(max(abs(u)), max(abs(u_0)))])


    %- check consistency of picks -----------------------------------------
    while true

        fprintf('\n')
        disp('pick left start of window:');
        [left, ~] = ginput(1);
        fprintf('left pick: %5.2f\n', left)
        fprintf('\n')
        disp('pick end of window:');
        [right, ~] = ginput(1);
        fprintf('right pick: %5.2f\n', right)

        if (strcmp(measurement, 'log_amplitude_ratio'))
            if (left * right < 0)
                warning('for log_amplitude_ratio measurement, pick on either causal or acausal branch')
                continue
            end
        end

        if (left < t(1) || left > t(end))
            warning('left pick is outside of bound')
            continue
        end

        if (right < t(1) || right > t(end))
            warning('right pick is outside of bound')
            continue
        end

        if (left >= right)
            warning('left pick is not left')
            continue
        end

        break

    end


    %- close figure -------------------------------------------------------
    close(fig)


end


