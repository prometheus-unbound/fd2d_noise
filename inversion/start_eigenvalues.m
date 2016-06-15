

addpath(genpath('../'))


usr_par.cluster = 'monch';
% 'local';
% 'monch';
% 'euler';
% 'brutus';


% usr_par.network = load( '../output/interferometry/array_1_ref_testing_small.mat' );
usr_par.network = load( '../output/interferometry/array_16_ref_hessian_small.mat' );


% set up options
opts.issym = 1;
opts.isreal = 1;
opts.disp = 2;


% start matlabpool
if( ~strcmp( usr_par.cluster, 'local') )
    parobj = start_cluster( usr_par.cluster, '', size(usr_par.network.ref_stat,1));
end


% get configuration and initialize counter
[~,~,nx,nz,~,~,~,~,~,n_basis_fct] = input_parameters();

global counter;
counter = 0;


% calculate eigenvalues
if( n_basis_fct == 0)
    [v, d, flag] = eigs( @eigenvalues, nx*nz, 10, 'lm', opts );
    % [v, d, flag] = eigs( @eigenvalues, nx*nz*2, 10, 'lm', opts );
else
    [v, d, flag] = eigs( @eigenvalues, nx*nz*n_basis_fct, 10, 'lm', opts );
    % [v, d, flag] = eigs( @eigenvalues, nx*nz*(n_basis_fct+1), 10, 'lm', opts );
end


% save solution
save( '../output/eigenvectors.mat', 'flag', 'v', 'd', 'counter' ,'-v7.3' )


% close matlabpool
if( ~strcmp( usr_par.cluster, 'local' ) )
    delete(parobj)
end