
addpath(genpath('../../'))


%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;

%% Define argument types for entry-point 'run_forward2_correlation'.
ARGS = cell(1,1);
ARGS{1} = cell(5,1);
ARGS{1}{1} = struct;
ARGS{1}{1}.rho = coder.typeof(0,[inf inf]);
ARGS{1}{1}.mu = coder.typeof(0,[inf inf]);
ARGS{1}{1} = coder.typeof(ARGS{1}{1});
ARGS{1}{2} = struct;
ARGS{1}{2}.spectrum = coder.typeof(0,[1 inf]);
ARGS{1}{2}.distribution = coder.typeof(0,[inf inf]);
ARGS{1}{2} = coder.typeof(ARGS{1}{2});
ARGS{1}{3} = coder.typeof(1i,[inf inf inf]);
ARGS{1}{4} = coder.typeof(0,[inf inf]);
ARGS{1}{5} = coder.typeof(0);

%% Invoke MATLAB Coder.
if( ~exist( 'exit_code', 'var' ) )
    exit_code = 0;
end

if( exit_code == 0 )
    try
        codegen -config cfg run2_forward_correlation -args ARGS{1}
        exit_code = 0;
    catch
        exit_code = 1;
    end
end

clear cfg
clear ARGS
rmpath(genpath('../../'))
