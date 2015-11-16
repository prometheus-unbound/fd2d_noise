
function [m] = map_parameters_to_m(m_parameters, usr_par)
% MAP_PARAMETERS_TO_M function to map the physical parameters to the model 
% variable m.
%
% Input: 
% m_parameters
% usr_par : auxiliary user defined parameters (optional)
%
% Output:
% m
%
% See also MAP_M_TO_PARAMETERS and MAP_GRADPARAMETERS_TO_GRADM.


if( strcmp( usr_par.type, 'source') )
    
    m = reshape( m_parameters, [], 1);
    
elseif( strcmp( usr_par.type, 'structure') )
    
%     m = reshape( sqrt( m_parameters ./ (usr_par.structure_inversion.v0^2 * usr_par.structure_inversion.rho ) ) - 1, [], 1 );
    m = reshape( m_parameters, [], 1);
    
end


end