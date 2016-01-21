
addpath(genpath('../../'))


% G_2, df
i1 = coder.typeof(complex(zeros(2,2,2)), [inf,inf,inf], 1);

% mu, ad-stf, spectrum
i2 = coder.typeof(zeros(2,2), [inf,inf], 1);

% src, rec
i3 = coder.typeof(zeros(2,2), [inf,2], 1);

% integer switch
i4 = coder.typeof(1, 1, 0);

% source_dist
i5 = coder.typeof(zeros(2,2,2), [inf,inf,inf], 1);

% forward_dxu_time, forward_dzu_time
i6 = coder.typeof((zeros(2,2,2,'single')), [inf,inf,inf], 1);


% run_forward1_green(mu, rho, src, mode)
codegen run_forward1_green.m -args {i2, i2, i3, i4}

% run_forward2_correlation(mu, rho, G_2, spectrum, source_dist, rec, mode, df)
codegen run_forward2_correlation.m -args {i2, i2, i1, i2, i5, i3, i4, i1}

% run_noise_source_kernel( mu, rho, G_2, spectrum, adstf, adsrc )
codegen run_noise_source_kernel.m -args {i2, i2, i1, i2, i2, i3}

% run_noise_structure_kernel( mu, rho, forward_dxu_time, forward_dzu_time, adstf, adsrc, spectrum, source_dist )
codegen run_noise_structure_kernel.m -args {i2, i2, i6, i6, i1, i3, i2, i5}


rmpath(genpath('../../'))