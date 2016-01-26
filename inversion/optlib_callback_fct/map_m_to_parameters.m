
function [m_parameters] = map_m_to_parameters(m, usr_par)
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


[~,~,nx,nz,~,~,~,~,~,n_basis_fct] = input_parameters();

if( n_basis_fct == 0 )
    n_basis_fct = 1;
end

m_parameters = imfilter( reshape( m, nx, nz, n_basis_fct + 1 ), usr_par.kernel.imfilter, 'circular' );
m_parameters(:,:,end) = 4.8e10 * ( 1 + m_parameters(:,:,end) );


end