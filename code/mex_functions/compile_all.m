
addpath(genpath('../../'))


% G_2, C_2_dxv, C_2_dzv, df
i1 = coder.typeof(complex(zeros(2,2,2)),[inf,inf,inf],1);

% mu, ad-stf, spectrum
i2 = coder.typeof(zeros(2,2),[inf,inf],1);

% src, rec
i3 = coder.typeof(zeros(2,2),[inf,2],1);

% integer switch
i4 = coder.typeof(1,1,0);

% source_dist
i5 = coder.typeof(zeros(2,2,2),[inf,inf,inf],1);


% run_forward_green_fast(mu,src)
codegen run_forward_green_fast.m -args {i2,i3}

% run_forward_correlation_fast(G_2, source_dist, noise_spectrum, mu, rec, mode, df)
codegen run_forward_correlation_fast.m -args {i1,i5,i2,i2,i3,i4,i1}

% run_noise_source_kernel_fast(G_2, mu,spectrum, stf, adsrc)
codegen run_noise_source_kernel_fast.m -args {i1,i2,i2,i2,i3}

% run_noise_mu_kernel_fast(C_2_dxv, C_2_dzv, mu, stf, adsrc)
codegen run_noise_mu_kernel_fast.m -args {i1,i1,i2,i2,i3}


rmpath(genpath('../../'))