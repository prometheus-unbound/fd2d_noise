
%==========================================================================
% [h] = plot_recordings( u, t, color, normalize )
%
% input:
%--------
% u: displacement recordings
% t: time axis
% color, e.g. 'r'
% normalize, true or false
%
% output:
%--------
% h: handle to plot of seismograms
%
%==========================================================================


function [h] = plot_recordings(u, t, color, normalize, array, ref_stat)


    spacing = 2;
    offset = 0;


    %- convert to velocity ------------------------------------------------
    % nt = length(t);
    % v = zeros( size(u,1), size(u,2), nt-1);
    %
    % for i_ref = 1:size(u,1)
    %     for i_rec = 1:size(u,2)
    %         v(i_ref,i_rec,:) = diff(u(i_ref,i_rec,:)) / (t(2)-t(1));
    %     end
    % end
    %
    % t = t(1:nt-1);
    % u = v;


    %- filter displacement recordings -------------------------------------
    for i = 1:size(u, 1)
        for j = 1:size(u, 2)
            u(i, j,:) = filter_seismogram(squeeze(u(i, j,:))', t, 0.02, 0.2, 4);
        end
    end
    
    
    %- sort by distance if array is given ---------------------------------
    if( exist('array','var') && exist('ref_stat','var') )
        
        dist = zeros( (size(array, 1)-1)* size(ref_stat, 1), 3 );
        for i = 1:size(ref_stat, 1)
            
            src = ref_stat(i,:);
            rec = array(~ismember(array, src, 'rows'),:);
            
            dist( 1+(i-1)*size(rec,1):i*size(rec,1), 1 ) = sqrt( (rec(:,1)-src(1,1)).^2 + (rec(:,2)-src(1,2)).^2 );
            dist( 1+(i-1)*size(rec,1):i*size(rec,1), 2 ) = i;
            dist( 1+(i-1)*size(rec,1):i*size(rec,1), 3 ) = 1:size(rec,1);
            
        end
        
        dist = sortrows( dist, 1 );
        
    else
        
        dist = zeros( size(u, 1) * size(u, 2), 3 );
        for i = 1:size(u, 1)
            dist(1+(i-1) * size(u, 2):i*size(u, 2),2) = i;
            dist(1+(i-1) * size(u, 2):i*size(u, 2),3) = 1:size(u, 2);
        end
        
    end


    %- plot recordings ----------------------------------------------------
    ax1 = gca;
    set(ax1, 'FontSize', 12)
    hold(ax1, 'on')
    k = 1;
    
    for i = 1:size(u, 1)

        for j = 1:size(u, 2)

            %- normalization ----------------------------------------------
            if (normalize == true)
                m = max(max(abs( u(dist(k,2), dist(k,3), : ) )));
                if( m == 0 )
                    m = 1;
                end
            else
                m = 1;
            end

            %- plot -------------------------------------------------------
            h = plot(ax1, t, spacing * k + squeeze( u(dist(k,2), dist(k,3), : ) ) / m + offset, ...
                'color', color, 'LineWidth', 1);
            k = k + 1;
            
        end

    end


    %- labels and limits --------------------------------------------------
    xlabel(ax1, 'time [s]');
    if (normalize)
        ylabel(ax1, 'normalized');
        % else
        %     ylabel('seismograms','FontSize',18);
    end

    ylimits = get(ax1, 'YLim');
    % set(ax1, 'YLim', [ylimits(1) - 0.2, ylimits(2) + 0.2])
    set(ax1, 'YLim', [ylimits(1) - 0.2, 0.2 + spacing * (k-1) + max(abs(squeeze(u(end, end,:))))/m] )
    set(ax1, 'YTick', []);


end


