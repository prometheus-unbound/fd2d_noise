
%==========================================================================
% [left,right] = pick_window( u, u_0, t, measurement, i )
%==========================================================================


function [left, right] = pick_window(u, u_0, t, measurement, i)


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
    ax1 = gca;
    set(ax1, 'FontSize', 12);
    hold(ax1, 'on')
    box(ax1, 'on')
    legend_string = [];


    %- plot seismograms ---------------------------------------------------
    plot(ax1, t, u, 'r')
    legend_string{end + 1} = 'synthetics';
    if( sum(u_0) ~= 0 )
        plot(ax1, t, u_0, 'k')
        legend_string{end + 1} = 'data';
    end


    % legend, labels, title, limits ---------------------------------------
    legend(ax1, legend_string)
    xlabel(ax1, 'time [s]')
    title(ax1, ['correlation ' num2str(i) ':   pick LEFT start of measurement window'], 'FontSize', 14)
    xlim(ax1, [t(1), t(end)])
    ylim(ax1, [- 1.2 * max(max(abs(u)), max(abs(u_0))), 1.2 * max(max(abs(u)), max(abs(u_0)))])


    %- check consistency of picks -----------------------------------------
    error = false;
    error_msg = [];
    while true

        fprintf('\n')
        disp('pick left start of window:');
        [left, dummy1] = ginput(1);
        h(1,:) = text(left,dummy1, sprintf('left pick:\n%5.2f s', left), 'FontSize',14);
        fprintf('left pick: %5.2f\n\n', left)
        
        title(ax1, ['correlation ' num2str(i) ':   pick END of measurement window'], 'FontSize', 14)
        disp('pick end of window:');
        [right, dummy2] = ginput(1);
        h(2,:) = text(right,dummy2, sprintf('right pick:\n%5.2f s', right),'FontSize',14);
        fprintf('right pick: %5.2f\n\n', right)

        if (strcmp(measurement, 'log_amplitude_ratio'))
            if (left * right < 0)
                error = true;
                error_msg{end + 1} = 'for log_amplitude_ratio measurement, pick on either causal or acausal branch\n';
            end
        end

        if (left < t(1) || left > t(end))
            error = true;
            error_msg{end + 1} = 'left pick is outside of bound\n';
        end

        if (right < t(1) || right > t(end))
            error = true;
            error_msg{end + 1} = 'right pick is outside of bound\n';
        end

        if (left >= right)
            error = true;
            error_msg{end + 1} = 'left pick is not left\n';
        end

        
        if( error )
            
            fprintf(['\n' error_msg{:}])
            
            dummy1 = get(gca,'XLim');
            dummy2 = get(gca,'YLim');
            h(3,:) = text( dummy1(1) + 0.05*abs(diff(dummy1)), dummy2(2) - 0.1*abs(diff(dummy2)), ...
                sprintf([error_msg{:}]), 'FontSize', 18, 'Color', 'r', 'Interpreter', 'none');
            
            pause(2.5)
            delete(h)
            error = false;
            error_msg = [];
            continue
            
        else
            break
        end

    end


    %- close figure -------------------------------------------------------
    pause(1)
    close(fig)


end


