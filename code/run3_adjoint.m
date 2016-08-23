
%==========================================================================
% compute sensitivity kernel for mu
%
% [K, u_adj_fft] = run3_adjoint( structure, noise_source, G_fft,
%       ref_station, rec, adjstf, wavefield_fwd, mode )
%
% input:
%--------
% structure: contains mu [N/m^2] and rho [kg/m^3]
% noise_source: contains spectrum and distribution of psd
% G_fft: Fourier transformed displacement Green function for ref. station
% ref_station: position of reference station
% rec: receiver positions, i.e. adjoint source locations
% adjstf: adjoint source time functions
% wavefield_fwd: forward wavefield
%
% mode: integer switch
%       == 0: do not Fourier transform adjoint wavefield
%       == 1: Fourier transform adjoint wavefield
%
% output:
%--------
% K: sensitivity kernels for source distribution and mu
% u_adj_fft: Fourier transformed adjoint wavefield - necessary for
%            structure kernel
%
%==========================================================================


function [K, u_adj_fft] = run3_adjoint(structure, noise_source, G_fft, ref_station, rec, adjstf, wavefield_fwd, mode)


    %- get basic configuration --------------------------------------------
    [Lx, Lz, nx, nz, dt, nt, order, ~, ~, store_fwd_nth, make_plots, plot_nth] = input_parameters();
    [X, Z, x, z, dx, dz] = define_computational_domain(Lx, Lz, nx, nz);


    %- time axis ----------------------------------------------------------
    t = - (nt - 1) * dt:dt:(nt - 1) * dt;
    n_zero = nt;
    nt = length(t);


    %- prepare coefficients for Fourier transform -------------------------
    [~, n_sample, w_sample, dw, freq_samp] = input_interferometry();

    fft_coeff = zeros(nt, n_sample) + 1i * zeros(nt, n_sample);
    ifft_coeff = zeros(nt, n_sample) + 1i * zeros(nt, n_sample);
    for n = 1:nt
        for k = 1:n_sample
            fft_coeff(n, k) = 1 / sqrt(2 * pi) * exp(- 1i * w_sample(k) * t(n)) * dt;
            ifft_coeff(n, k) = 1 / sqrt(2 * pi) * exp(1i * w_sample(k) * t(n)) * dw;
        end
    end


    %- compute indices for receivers, i.e. the adjoint source locations ---
    n_receivers = size(rec, 1);
    rec_id = zeros(n_receivers, 2);

    for i = 1:n_receivers
        rec_id(i, 1) = find( min( abs(x - rec(i, 1)) ) == abs(x - rec(i, 1)), 1 );
        rec_id(i, 2) = find( min( abs(z - rec(i, 2)) ) == abs(z - rec(i, 2)), 1 );
    end


    %- initialise absorbing boundary taper a la Cerjan --------------------
    [absbound] = init_absbound();


    %- prepare figure for kernel building process -------------------------
    if (strcmp(make_plots, 'yes') && isempty(wavefield_fwd))

        fig = figure;
        set(fig, 'units', 'normalized', 'position', [0.1, 0.3, 0.6, 0.5])

        ax1 = subplot(1, 2, 1);
        hold on
        set(ax1, 'FontSize', 18);
        xlabel('x [km]')
        ylabel('z [km]')
        set(ax1, 'XTick', [0, 200, 400])
        set(ax1, 'YTick', [0, 200, 400])
        title('forward and adjoint wavefield', 'FontSize', 22)
        cm = cbrewer('div', 'RdBu', 120, 'PCHIP');
        colormap(cm)
        axis square
        box on
        set(gca, 'LineWidth', 2)


        ax2 = subplot(1, 2, 2);
        hold on
        set(ax2, 'FontSize', 18);
        xlabel('x [km]')
        set(ax2, 'XTick', [0, 200, 400])
        set(ax2, 'YTick', [])
        title('kernel build-up', 'FontSize', 22)
        colormap(cm)
        axis square
        box on
        set(gca, 'LineWidth', 2)

        cb1 = colorbar('Position', [0.50, 0.34, 0.02, 0.37], ...
            'Ticks', [], 'AxisLocation', 'in');
        ylabel(cb1, 'fields and kernels are normalized')
        cb2 = colorbar('Position', [0.50, 0.34, 0.02, 0.37], ...
            'Ticks', [- 1, 1], 'TickLabels', {'-', '+'}, 'AxisLocation', 'out');

        [width] = absorb_specs();
        max_u = 0;
        max_M_tn = 0;

    end


    %- allocate fields ----------------------------------------------------
    v = zeros(nx, nz);
    sxy = zeros(nx - 1, nz);
    szy = zeros(nx, nz - 1);
    u = zeros(nx, nz);

    K_mu = zeros(nx, nz);
    K_source = zeros(nx, nz);

    i_fw_in = 1;
    if (mode == 1)
        u_adj_fft = zeros(nx, nz, n_sample) + 1i * zeros(nx, nz, n_sample);
    else
        u_adj_fft = [];
    end


    %======================================================================
    % time loop
    %======================================================================

    %- only need first half of second adjoint run
    if (size(adjstf, 3) ~= 1)
        nt = n_zero;
    end


    for n = 1:nt


        %- compute divergence of current stress tensor --------------------
        DS = div_s(sxy, szy, dx, dz, nx, nz, order);


        %- add adjoint source time function -------------------------------
        if (size(adjstf, 3) == 1 && ~isempty(adjstf))

            for i = 1:n_receivers
                DS(rec_id(i, 1), rec_id(i, 2)) = DS(rec_id(i, 1), rec_id(i, 2)) + real(adjstf(i, n));
            end

        elseif (size(adjstf, 3) ~= 1 && ~isempty(adjstf))

            if (mod(n, freq_samp) == 0)
                T = zeros(nx, nz) + 1i * zeros(nx, nz);

                for k = 1:n_sample
                    T = T + noise_source.spectrum(k) * noise_source.distribution ...
                        .* conj(adjstf(:,:,k)) * ifft_coeff(n, k);
                end

                DS = DS + real(T);

            end

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


        %- calculate displacement and respective strain -------------------
        u = u + v * dt;
        strain_dxu = dx_v(u, dx, nx, nz, order);
        strain_dzu = dz_v(u, dz, nx, nz, order);


        %- build up source kernel -----------------------------------------
        if (~isempty(G_fft) && size(G_fft, 3) == n_sample && mod(n, freq_samp) == 0)

            M_tn = zeros(nx, nz) + 1i * zeros(nx, nz);

            for k = 1:n_sample
                M_tn = M_tn + noise_source.spectrum(k) .* G_fft(:,:,k) * ifft_coeff(n, k);
            end

            K_source = K_source + real(M_tn .* u);

        end


        %- build up structure kernel --------------------------------------
        if (~isempty(wavefield_fwd) && size(wavefield_fwd, 3) >= i_fw_in && mod(n, store_fwd_nth) == 0)

            K_mu(1:nx - 1,:) = K_mu(1:nx-1, :) - ...
                strain_dxu .* dx_v(wavefield_fwd(:,:, end - i_fw_in + 1), dx, nx, nz, order) * store_fwd_nth;

            K_mu(:, 1:nz - 1) = K_mu(:, 1:nz-1) - ...
                strain_dzu .* dz_v(wavefield_fwd(:,:, end - i_fw_in + 1), dz, nx, nz, order) * store_fwd_nth;

            i_fw_in = i_fw_in + 1;

        end


        %- save Fourier transformed adjoint state for second run ----------
        if (mode == 1)

            if (mod(n, freq_samp) == 0)
                for k = 1:n_sample
                    u_adj_fft(:,:, k) = u_adj_fft(:,:, k) + u * fft_coeff(n, k);
                end
            end

        end


        %- plot correlation wavefield -------------------------------------
        if (strcmp(make_plots, 'yes') && isempty(wavefield_fwd))

            if (mod(n, plot_nth) == 0 || n == nt)

                subplot(1, 2, 1)
                cla

                max_u = max(max_u, max(max(abs(u + eps))));
                max_M_tn = max(max_M_tn, max(max(abs(real(M_tn + eps)))));

                if (t(n) >= 0)
                    pcolor(X / 1000, Z / 1000, u' / max_u + real(M_tn)' / max_M_tn);
                elseif any(max(adjstf(:, 1:n), [], 2) > 0.05 * max(adjstf(:, 1:n_zero), [], 2))
                    pcolor(X / 1000, Z / 1000, u' / max_u);
                else
                    pcolor(X / 1000, Z / 1000, 0 * u');
                end
                caxis([- 0.2, 0.2])

                plot(ref_station(:, 1) / 1000, ref_station(:, 2) / 1000, ...
                    'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
                plot(rec(:, 1) / 1000, rec(:, 2) / 1000, ...
                    'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8)

                plot([width, Lx - width] / 1000, [width, width] / 1000, 'k--');
                plot([width, Lx - width] / 1000, [Lz - width, Lz - width] / 1000, 'k--')
                plot([width, width] / 1000, [width, Lz - width] / 1000, 'k--')
                plot([Lx - width, Lx - width] / 1000, [width, Lz - width] / 1000, 'k--')

                shading interp


                subplot(1, 2, 2)
                if (t(n) >= 0)
                    pcolor(X / 1000, Z / 1000, K_source')
                    m = max(max(abs(K_source)));
                    caxis([- 0.6 * m, 0.6 * m]);
                else
                    pcolor(X / 1000, Z / 1000, 0 * K_source');
                end

                plot(ref_station(:, 1) / 1000, ref_station(:, 2) / 1000, ...
                    'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
                plot(rec(:, 1) / 1000, rec(:, 2) / 1000, ...
                    'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8)

                plot([width, Lx - width] / 1000, [width, width] / 1000, 'k--');
                plot([width, Lx - width] / 1000, [Lz - width, Lz - width] / 1000, 'k--')
                plot([width, width] / 1000, [width, Lz - width] / 1000, 'k--')
                plot([Lx - width, Lx - width] / 1000, [width, Lz - width] / 1000, 'k--')
                xlim([0, Lx / 1000])
                ylim([0, Lz / 1000])

                set(cb2, 'Ticks', get(cb2, 'Limits'))

                shading interp
                drawnow

            end

        end


    end


    %- concatenate kernels ------------------------------------------------
    K = zeros(nx, nz, 2);
    K(:,:,1) = K_mu;
    K(:,:,2) = K_source;


end

