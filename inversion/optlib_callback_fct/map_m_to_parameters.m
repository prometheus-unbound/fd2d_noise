
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


m_parameters = reshape( m, usr_par.config.nx, usr_par.config.nz, [] );

[absbound] = init_absbound();
[ix,iz] = find( absbound == 1, 1, 'first' );

m_parameters(ix:end-ix+1, iz:end-iz+1, 1:end-1) = imfilter( m_parameters(ix:end-ix+1, iz:end-iz+1, 1:end-1) , usr_par.kernel.imfilter.source, 'circular' );
m_parameters(ix:end-ix+1, iz:end-iz+1, end) = imfilter( m_parameters(ix:end-ix+1, iz:end-iz+1, end) , usr_par.kernel.imfilter.structure, 'circular' );

m_parameters(:,:,1:end-1) = usr_par.initial.ref_source + m_parameters(:,:,1:end-1);
m_parameters(:,:,end) = usr_par.initial.mu_0 * ( usr_par.initial.ref_structure + m_parameters(:,:,end) );


end