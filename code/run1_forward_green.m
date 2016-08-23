
%==========================================================================
% compute Green function for reference station
%
% [ G_fft, G_out ] = run1_forward_green( structure, ref_station, mode )
%
% input:
%--------
% structure: contains mu [N/m^2] and rho [kg/m^3]
% ref_station: position of reference station
% mode: integer switch
%       == 0: do not save wavefield
%       == 1: save wavefield
%
% output:
%--------
% G_fft: Fourier transformed displacement Green function for ref. station
% G_out: wavefield of displacement Green function
%
%==========================================================================


function [G_fft, G_out] = run1_forward_green(structure, ref_station, mode)


    %- get basic configuration --------------------------------------------
    [Lx, Lz, nx, nz, dt, nt, order, ~, ~, store_fwd_nth] = input_parameters();
    [~, ~, x, z, dx, dz] = define_computational_domain(Lx, Lz, nx, nz);


    %- time axis ----------------------------------------------------------
    t = 0:dt:(nt - 1) * dt;


    %- prepare coefficients for Fourier transform -------------------------
    [~, n_sample, w_sample, ~, freq_samp] = input_interferometry();

    fft_coeff = zeros(nt, n_sample) + 1i * zeros(nt, n_sample);
    for n = 1:nt
        for k = 1:n_sample
            fft_coeff(n, k) = 1 / sqrt(2 * pi) * exp(- 1i * w_sample(k) * t(n)) * dt;
        end
    end


    %- make source time function ------------------------------------------
    stf = 1.0e9 * ones(1, nt);


    %- compute indices for source locations -------------------------------
    assert(size(ref_station, 1) == 1, ...
        'not possible to compute Green functions for more reference stations in one run of run1_forward_green.m')
    src_id = zeros(1, 2);
    src_id(1, 1) = find( min( abs(x - ref_station(1, 1)) ) == abs(x - ref_station(1, 1)), 1 );
    src_id(1, 2) = find( min( abs(z - ref_station(1, 2)) ) == abs(z - ref_station(1, 2)), 1 );


    %- initialise absorbing boundary taper a la Cerjan --------------------
    [absbound] = init_absbound();


    %- allocate fields ----------------------------------------------------
    v = zeros(nx, nz);
    sxy = zeros(nx - 1, nz);
    szy = zeros(nx, nz - 1);
    G_fft = zeros(nx, nz, n_sample) + 1i * zeros(nx, nz, n_sample);

    i_fwd_out = 1;
    n_fwd = 0;
    for n = nt:(2 * nt - 1)
        if (mod(n, store_fwd_nth) == 0)
            n_fwd = n_fwd + 1;
        end
    end

    if (mode ~= 0)
        G_out = zeros(nx, nz, n_fwd, 'single');
        % G_out = zeros(nx,nz,n_fw);
    else
        G_out = single([]);
    end


    %======================================================================
    % time loop
    %======================================================================

    for n = 1:nt


        if (mode ~= 0 && mod(n + nt - 1, store_fwd_nth) == 0)
            G_out(:,:, i_fwd_out) = single(v);
            % G_out(:,:,i_fw_out) = v;
            i_fwd_out = i_fwd_out + 1;
        end


        %- compute divergence of current stress tensor --------------------
        DS = div_s(sxy, szy, dx, dz, nx, nz, order);


        %- add point source -----------------------------------------------
        DS(src_id(1, 1), src_id(1, 2)) = DS(src_id(1, 1), src_id(1, 2)) + stf(n);


        %- update velocity field ------------------------------------------
        v = v + dt * DS ./ structure.rho;


        %- apply absorbing boundary taper ---------------------------------
        v = v .* absbound;


        %- compute derivatives of current velocity and update stress tensor
        strain_dxv = dx_v(v, dx, nx, nz, order);
        strain_dzv = dz_v(v, dz, nx, nz, order);

        sxy = sxy + dt * structure.mu(1:nx-1, :) .* strain_dxv;
        szy = szy + dt * structure.mu(:, 1:nz-1) .* strain_dzv;


        %- accumulate Fourier transform of the displacement Greens function
        if (mod(n + nt - 1, freq_samp) == 0)

            for k = 1:n_sample
                G_fft(:,:,k) = G_fft(:,:,k) + v * fft_coeff(n, k);
            end

        end


    end


end


