
function [m] = map_parameters_to_m( m_parameters, usr_par )
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


m = reshape( m_parameters, [], 1 );
    

end