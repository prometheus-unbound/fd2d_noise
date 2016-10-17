
addpath(genpath('../../'))


%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;

%% Define argument types for entry-point 'run_forward1_green'.
ARGS = cell(1, 1);
ARGS{1} = cell(3, 1);
ARGS{1}{1} = struct;
ARGS{1}{1}.rho = coder.typeof(0, [inf, inf]);
ARGS{1}{1}.mu = coder.typeof(0, [inf, inf]);
ARGS{1}{1} = coder.typeof(ARGS{1}{1});
ARGS{1}{2} = coder.typeof(0, [inf, 2]);
ARGS{1}{3} = coder.typeof(0);

%% Invoke MATLAB Coder.
try
    codegen -config cfg run1_forward_green -args ARGS{1}
    exit_code = 0;
catch
    exit_code = 1;
end

clear cfg
clear ARGS
rmpath(genpath('../../'))

