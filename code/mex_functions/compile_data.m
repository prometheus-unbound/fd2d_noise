
addpath(genpath('../../'))


% G_2
i1 = coder.typeof(complex(zeros(2,2,2)), [inf,inf,inf], 1);

% mu, rho, spectrum
i2 = coder.typeof(zeros(2,2), [inf,inf], 1);

% src, rec
i3 = coder.typeof(zeros(2,2), [inf,2], 1);

% integer switch, i.e. mode
i4 = coder.typeof(1, 1, 0);

% source_dist
i5 = coder.typeof(zeros(2,2,2), [inf,inf,inf], 1);



% run_forward1_green(mu, rho, src, mode)
codegen run_forward1_green.m -args {i2, i2, i3, i4}

% run_forward2_correlation(mu, rho, G_2, spectrum, source_dist, rec, mode)
codegen run_forward2_correlation.m -args {i2, i2, i1, i2, i5, i3, i4}



rmpath(genpath('../../'))