
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
        set(fig, 'units', 'normalized', 'position', [0.1, 0.3, 0.3, 0.5])
        hold on

        xlabel('x [km]')
        ylabel('z [km]')
        set(gca, 'XTick', [0, 200, 400])
        set(gca, 'YTick', [0, 200, 400])

        title('correlation wavefield', 'FontSize', 14)
        cm = cbrewer('div', 'RdBu', 120);
        colormap(cm)
        cb = colorbar;
        axis square
        box on
        set(gca, 'LineWidth', 2)

        [width] = absorb_specs();

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
        % if( mod(n,freq_samp) == 0 && t(n) <= 0.0 )
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

                pcolor(X / 1000, Z / 1000, u')
                shading interp

                if (strcmp(source_type, 'homogeneous'))
                    caxis([-0.08, 0.08]);
                elseif (strcmp(source_type, 'point'))
                    caxis([-0.003, 0.003]);
                elseif (strcmp(source_type, 'gaussian'))
                    caxis([-0.2, 0.2]);
                else
                    m = max(max(abs(u)));
                    caxis([- 0.8 * m, 0.8 * m]);
                end

                clabels = get(cb, 'YTick');
                set(cb, 'YTick', [clabels(1), clabels(ceil(length(clabels) / 2)), clabels(end)])

                plot(ref_station(:, 1) / 1000, ref_station(:, 2) / 1000, ... 
                    'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
                plot(rec(:, 1) / 1000, rec(:, 2) / 1000, ....
                    'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8)

                plot([width, Lx - width] / 1000, [width, width] / 1000, 'k--');
                plot([width, Lx - width] / 1000, [Lz - width, Lz - width] / 1000, 'k--')
                plot([width, width] / 1000, [width, Lz - width] / 1000, 'k--')
                plot([Lx - width, Lx - width] / 1000, [width, Lz - width] / 1000, 'k--')
                
                set(gca, 'FontSize', 12, 'position', [0.17, 0.204, 0.599, 0.624]);
                drawnow

            end

        end


    end


    if (strcmp(make_plots, 'yes'))
        close(fig)
    end


end


