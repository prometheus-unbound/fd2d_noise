
%==========================================================================
% compute adjoint sources
%
% [misfit, adjstf] = make_adjoint_sources( u, u_0, t, measurement, src, rec, mode )
%
% input:
%-------
% u: synthetic displacement seismograms
% u_0: "observed" displacement seismograms
% t: time axis
% measurement:  'log_amplitude_ratio'
%               'amplitude_difference'
%               'waveform_difference'
%               'cc_time_shift'
% src: source position
% rec: receiver positions
%
% output:
%-------
% misfit: vector containing misfit of each receiver
% adstf: vector with adjoint source time function for each receiver
%
%
% When u_0, i.e. the observed displacement seismograms, are zero,
% the code performs data-independent measurements.
%
%==========================================================================


function [misfit, adjstf] = make_adjoint_sources(u, u_0, t, measurement, src, rec, mode)


    %- initialisations ----------------------------------------------------
    nt = length(t);
    n_rec = size(rec, 1);
    misfit = zeros(n_rec, 1);
    adjstf = zeros(n_rec, nt);


    %======================================================================
    %- march through the various recordings
    %======================================================================

    for i_rec = 1:n_rec


        %- select time windows --------------------------------------------
        if( strcmp(mode, 'manual') )

            [left, right] = pick_window(u(i_rec,:), u_0(i_rec,:), t, measurement, i_rec);

        else

            if (strcmp(measurement, 'waveform_difference'))
                left = t(1);
                right = t(end);

            else

                distance = sqrt((src(1, 1) - rec(i_rec, 1)) .^ 2 + (src(1, 2) - rec(i_rec, 2)) .^ 2);
                left = distance / 4000.0 - 27.0;
                right = distance / 4000.0 + 27.0;

                if (left < 0)
                    index = find(t == 0);
                    left = t(index + 1);
                end

                if (right > t(end))
                    right = t(end);
                end

                %- error if time series is not long enough ----------------
                if (left > t(end))
                    error('time series too short?')
                end

            end

        end


        %- define window function -----------------------------------------
        win = get_window(t, left, right, 'cos_taper');


        %- compute misfit and adjoint source time function ----------------
        if (strcmp(measurement, 'waveform_difference'))

            [misfit(i_rec,:), adjstf(i_rec,:)] = waveform_difference(u(i_rec,:), u_0(i_rec,:), win, t);


        elseif (strcmp(measurement, 'cc_time_shift'))

            [misfit_irec_caus, adjstf_irec_caus(1,:)] = cc_time_shift(u(i_rec,:), u_0(i_rec,:), win, t);

            [left, right] = swap(- left, - right);
            win = get_window(t, left, right, 'cos_taper');
            [misfit_irec_acaus, adjstf_irec_acaus(1,:)] = cc_time_shift(u(i_rec,:), u_0(i_rec,:), win, t);

            misfit(i_rec,:) = misfit_irec_caus + misfit_irec_acaus;
            adjstf(i_rec,:) = adjstf_irec_caus + adjstf_irec_acaus;


        elseif (strcmp(measurement, 'amplitude_difference'))

            % [misfit_irec_caus, adjstf_irec_caus(1,:)] = amp_diff( u(i_rec,:), u_0(i_rec,:), win, t );
            %
            % [left, right] = swap( -left, -right );
            % win = get_window( t, left, right, 'cos_taper' );
            % [misfit_irec_acaus, adjstf_irec_acaus(1,:)] = amp_diff( u(i_rec,:), u_0(i_rec,:), win, t );

            % misfit(i_rec,:) = misfit_irec_caus + misfit_irec_acaus;
            % adjstf(i_rec,:) = adjstf_irec_caus + adjstf_irec_acaus;

            % [left, right] = swap( -left, -right );
            % win = get_window( t, left, right, 'cos_taper' );
            [misfit(i_rec,:), adjstf(i_rec,:)] = amp_diff(u(i_rec,:), u_0(i_rec,:), win, t);


        elseif (strcmp(measurement, 'log_amplitude_ratio'))

            % win = get_window(t, left, right, 'hann');
            [misfit(i_rec,:), adjstf(i_rec,:)] = log_amp_ratio(u(i_rec,:), u_0(i_rec,:), win, t);


        end


    end


end


