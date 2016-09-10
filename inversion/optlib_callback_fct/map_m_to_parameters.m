
function [m_parameters] = map_m_to_parameters( m, usr_par )
% MAP_M_TO_PARAMETERS function to map the model variable m to physical
% parameters.
%
% Input: 
% m
% usr_par : auxiliary user defined parameters (optional)
%
% Output:
% m_parameters
%
% See also MAP_PARAMETERS_TO_M and MAP_GRADPARAMETERS_TO_GRADM.


[Lx, Lz, nx, nz] = input_parameters();
[~, ~, x, z] = define_computational_domain(Lx, Lz, nx, nz);

m_parameters = reshape( m, nx, nz, [] );

[absbound] = init_absbound();
[ix,iz] = find( absbound == 1, 1, 'first' );

m_parameters(ix:end-ix+1, iz:end-iz+1, 1) = gaussblur2d(m_parameters(ix:end-ix+1, iz:end-iz+1, 1), x(ix:end-ix+1), z(iz:end-iz+1), usr_par.smoothing.sigma);
m_parameters(ix:end-ix+1, iz:end-iz+1, 2) = gaussblur2d(m_parameters(ix:end-ix+1, iz:end-iz+1, 2), x(ix:end-ix+1), z(iz:end-iz+1), usr_par.smoothing.sigma);


end