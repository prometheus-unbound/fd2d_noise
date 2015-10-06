
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


[~,~,nx,nz] = input_parameters();


if( strcmp( usr_par.type, 'source') )
    
    m_parameters = reshape( m, nx, nz );
    
elseif( strcmp( usr_par.type, 'structure') )
    
    % in this case, m_parameters is mu, we don't consider rho at the moment
    m_parameters = reshape( usr_par.structure_inversion.v0^2 * reshape( usr_par.structure_inversion.rho, [], 1) .* (1+m).^2, nx, nz );
    
    % m_parameters = reshape( m, nx, nz );
    
end


end