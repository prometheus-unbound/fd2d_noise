
%==========================================================================
% generate noise source distribution
%
% [ noise_source ] = make_noise_source( make_plots )
%
% input:
%--------
% make_plots: 'yes' or 'no' (optional)
%
% output:
%--------
% noise_source: contains spectrum and distribution of psd
%
%==========================================================================


function [noise_source] = make_noise_source(make_plots)


    %======================================================================
    % user input
    %======================================================================

    % specify spectrum
    f_peak = 1 / 9;
    bandwidth = 0.04;
    strength = 1;

    % specify point location, width and strength
    point.x_source = 2.0e5;
    point.z_source = 1.7e5;
    point.source_width = 5e3;
    point.magnitude = 5.0;

    % specify location, width and strength of Gaussian anomaly
    gaussian.x_source = 1.0e5;
    gaussian.z_source = 1.6e5;
    gaussian.source_width = 3e4;
    gaussian.magnitude = 15.0;


    %======================================================================
    % set up noise source as specified by user
    %======================================================================

    %- define make_plots if not specified ---------------------------------
    if (nargin < 1)
        make_plots = 'yes';
    end


    %- get basic configuration --------------------------------------------
    [Lx, Lz, nx, nz, ~, ~, ~, ~, source_type] = input_parameters();
    [X, Z] = define_computational_domain(Lx, Lz, nx, nz);
    f_sample = input_interferometry();


    %- define source spectrum ---------------------------------------------
    noise_source.spectrum = strength * exp(- (abs(f_sample) - f_peak) .^ 2 / bandwidth ^ 2);


    %- define geographic power-spectral density distribution --------------
    if (strcmp(source_type, 'homogeneous'))

        noise_source.distribution = ones(nx, nz);

    elseif (strcmp(source_type, 'point'))

        noise_source.distribution = point.magnitude * ...
            exp(- ((X - point.x_source) .^ 2 + (Z - point.z_source) .^ 2) / point.source_width ^ 2)';

    elseif (strcmp(source_type, 'gaussian'))

        noise_source.distribution = ones(nx, nz);
        noise_source.distribution = noise_source.distribution + ...
            gaussian.magnitude * exp(- ((X - gaussian.x_source) .^ 2 + (Z - gaussian.z_source) .^ 2) / gaussian.source_width ^ 2)';

    end


    %- plot noise source configuration ------------------------------------
    if (strcmp(make_plots, 'yes'))

        figure
        ax1 = gca;
        set(ax1, 'FontSize', 12)
        hold(ax1, 'on')
        grid(ax1, 'on')
        plot(ax1, f_sample, noise_source.spectrum, 'r')
        xlabel(ax1, 'frequency [Hz]');
        xlim(ax1, [f_sample(1), f_sample(end)])
        title(ax1, 'spectrum for noise source', 'FontSize', 14)

        if (exist([fd2d_path(), 'output', filesep, 'array_1_ref.mat'], 'file'))
            load([fd2d_path(), 'output', filesep, 'array_1_ref.mat']);
        else
            array = [];
        end

        plot_models([], noise_source.distribution, array, [0, 0, 0, 0]);

    end


end


