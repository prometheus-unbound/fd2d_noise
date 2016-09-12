
%==========================================================================
% compute correlation wavefield
%
% [ seismograms, C_out ] = run2_forward_correlation( ...
%       structure, noise_source, G_fft, ref_station, rec, mode )
%
% input:
%--------
% structure: contains mu [N/m^2] and rho [kg/m^3]
% noise_source: contains spectrum and distribution of psd
% G_fft: Fourier transformed displacement Green function for ref. station
% ref_station: position of reference station
% rec: receiver positions
% mode: integer switch
%       == 0: do not save wavefield
%       == 1: save wavefield
%
% output:
%--------
% seismograms, i.e. correlation recordings
% C: correlation wavefield
%
%==========================================================================


function [seismograms, C_out] = run2_forward_correlation(structure, noise_source, G_fft, ref_station, rec, mode)


    %- get basic configuration --------------------------------------------
    [Lx, Lz, nx, nz, dt, nt, order, ~, source_type, store_fwd_nth, make_plots, plot_nth] = input_parameters();
    [X, Z, x, z, dx, dz] = define_computational_domain(Lx, Lz, nx, nz);


    %- time axis ----------------------------------------------------------
    t = - (nt - 1) * dt:dt:(nt - 1) * dt;
    nt = length(t);


    %- prepare coefficients for Fourier transform -------------------------
    [~, n_sample, w_sample, dw, freq_samp] = input_interferometry();

    ifft_coeff = zeros(nt, n_sample) + 1i * zeros(nt, n_sample);
    for n = 1:nt
        for k = 1:n_sample
            ifft_coeff(n, k) = 1 / sqrt(2 * pi) * exp(1i * w_sample(k) * t(n)) * dw;
        end
    end


    %- compute indices for receiver locations -----------------------------
    n_receivers = size(rec, 1);
    rec_id = zeros(n_receivers, 2);

    for i = 1:n_receivers
        rec_id(i, 1) = find( min( abs(x - rec(i, 1)) ) == abs(x - rec(i, 1)), 1 );
        rec_id(i, 2) = find( min( abs(z - rec(i, 2)) ) == abs(z - rec(i, 2)), 1 );
    end


    %- initialise absorbing boundary taper a la Cerjan --------------------
    [absbound] = init_absbound();


    %- prepare figure for correlation wavefield ---------------------------
    if (strcmp(make_plots, 'yes'))

        fig = figure;
        set(fig, 'units', 'normalized', 'position', [0.1, 0.3, 0.3, 0.5], 'Color', [1 1 1])
        title_size = 14;
        font_size = 12;
        marker_size = 6;
        ax1 = gca;
        hold(ax1, 'on')

        xlabel(ax1, 'x [km]')
        ylabel(ax1, 'z [km]')

        title(ax1, 'correlation wavefield', 'FontSize', title_size)
        cm = cbrewer('div', 'RdBu', 120);
        colormap(cm)
        cb = colorbar;
        axis(ax1, 'image')
        box(ax1, 'on')
        set(ax1, 'LineWidth', 2)

        [width, absorb_left, absorb_right, absorb_top, absorb_bottom] = absorb_specs();
        
        %- make movie -----------------------------------------------------
        % writerObj = VideoWriter('~/Desktop/test','MPEG-4');
        % writerObj.FrameRate = 6;
        % open(writerObj);

    end


    %- allocate fields ----------------------------------------------------
    v = zeros(nx, nz);
    sxy = zeros(nx - 1, nz);
    szy = zeros(nx, nz - 1);
    u = zeros(nx, nz);
    seismograms = zeros(n_receivers, nt);

    i_fwd_out = 1;
    n_fwd = floor(nt / store_fwd_nth);
    if (mode ~= 0)
        C_out = zeros(nx, nz, n_fwd, 'single');
    else
        C_out = single([]);
    end


    %======================================================================
    % time loop
    %======================================================================

    for n = 1:nt


        %- save correlation wavefield -------------------------------------
        if (mode ~= 0 && mod(n, store_fwd_nth) == 0)
            C_out(:,:, i_fwd_out) = single(u);
            i_fwd_out = i_fwd_out + 1;
        end


        %- compute divergence of current stress tensor --------------------
        DS = div_s(sxy, szy, dx, dz, nx, nz, order);


        %- add source of the correlation field ----------------------------
        if (mod(n, freq_samp) == 0)

            %- transform on the fly to the time domain
            S = zeros(nx, nz) + 1i * zeros(nx, nz);

            % calculate source for correlation wavefield
            for k = 1:n_sample
                S = S + noise_source.spectrum(k) * noise_source.distribution .* conj(G_fft(:,:,k)) * ifft_coeff(n, k);
            end

            DS = DS + real(S);

        end


        %- update velocity field ------------------------------------------
        v = v + dt * DS ./ structure.rho;


        %- apply absorbing boundary taper ---------------------------------
        v = v .* absbound;


        %- compute derivatives of current velocity and update stress tensor
        strain_dxv = dx_v(v, dx, nx, nz, order);
        strain_dzv = dz_v(v, dz, nx, nz, order);

        sxy = sxy + dt * structure.mu(1:nx-1, :) .* strain_dxv;
        szy = szy + dt * structure.mu(:, 1:nz-1) .* strain_dzv;


        %- calculate displacement -----------------------------------------
        u = u + v * dt;


        %- record seismograms ---------------------------------------------
        for ir = 1:n_receivers
            seismograms(ir, n) = u(rec_id(ir, 1), rec_id(ir, 2));
        end


        %- plot correlation wavefield -------------------------------------
        if (strcmp(make_plots, 'yes'))

            if (n > 0.1 * nt && n < 0.9 * nt && mod(n, plot_nth) == 0)

                
                %- plot wavefield -----------------------------------------
                pcolor(ax1, X / 1000, Z / 1000, u')
                shading(ax1, 'interp')

                if (strcmp(source_type, 'homogeneous'))
                    m = 0.08;
                elseif (strcmp(source_type, 'point'))
                    m = 0.003;                    
                elseif (strcmp(source_type, 'gaussian'))
                    m = 0.15;
                else
                    m = 0.8*max(max(abs(u)));
                end
                caxis(ax1, [-m, m]);

                
                %- plot array ---------------------------------------------
                plot(ax1, ref_station(:, 1) / 1000, ref_station(:, 2) / 1000, ... 
                    'kx', 'MarkerFaceColor', 'k', 'MarkerSize', marker_size)
                plot(ax1,rec(:, 1) / 1000, rec(:, 2) / 1000, ....
                    'kd', 'MarkerFaceColor', 'k', 'MarkerSize', marker_size)
                
                
                %- plot absorbing boundaries ------------------------------
                if(absorb_bottom); plot(ax1, [absorb_left*width, Lx - absorb_right*width] / 1000, [width, width] / 1000, 'k--'); end
                if(absorb_top); plot(ax1, [absorb_left*width, Lx - absorb_right*width] / 1000, [Lz - width, Lz - width] / 1000, 'k--'); end
                if(absorb_left); plot(ax1, [width, width] / 1000, [absorb_bottom*width, Lz - absorb_top*width] / 1000, 'k--'); end
                if(absorb_right); plot(ax1, [Lx - width, Lx - width] / 1000, [absorb_bottom*width, Lz - absorb_top*width] / 1000, 'k--'); end
                
                
                %- set labels ---------------------------------------------
                if( exist('OCTAVE_VERSION', 'builtin' ) == 0 )
                    set(cb, 'YTick', [-m m], 'TickLabels', {'-', '+'})
                else
                    clabels = get(cb, 'YTick');
                    set(cb, 'YTick', [clabels(1), clabels(ceil(length(clabels) / 2)), clabels(end)])
                end
                
                xlabels = get(ax1, 'XTick');
                ylabels = get(ax1, 'YTick');
                set(ax1, 'XTick', [xlabels(1), xlabels(ceil(length(xlabels) / 2)), xlabels(end)])
                set(ax1, 'YTick', [ylabels(1), ylabels(ceil(length(ylabels) / 2)), ylabels(end)])
                xlim(ax1, [0, Lx / 1000])
                ylim(ax1, [0, Lz / 1000])
                
                
                % set FontSize and invoke plot ----------------------------
                set(ax1, 'FontSize', font_size, 'position', [0.17, 0.204, 0.599, 0.624]);
                drawnow
                
                
                %- make movie ---------------------------------------------
                % M = getframe(fig);
                % writeVideo(writerObj,M);
                

            end

        end


    end

    
    if (strcmp(make_plots, 'yes'))
        % close(writerObj);
        close(fig)
    end


end


