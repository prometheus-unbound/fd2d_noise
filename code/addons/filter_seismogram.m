
%==========================================================================
% [u] = filter_seismogram(u, t, f_min, f_max, n_times)
%==========================================================================


function [u] = filter_seismogram(u, t, f_min, f_max, n_times)


    for j = 1:n_times

        u = fliplr(butterworth_lp(fliplr(u), t, 5, f_max, 'silent'));
        u = butterworth_lp(u, t, 5, f_max, 'silent');

        u = fliplr(butterworth_hp(fliplr(u), t, 3, f_min, 'silent'));
        u = butterworth_hp(u, t, 3, f_min, 'silent');

    end
    

end
