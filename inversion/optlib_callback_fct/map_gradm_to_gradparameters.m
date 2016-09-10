
function [grad_parameters] = map_gradm_to_gradparameters( m, grad_m, usr_par )
% MAP_GRADM_TO_GRADPARAMETERS function to map the model gradient gm to 
% physical parameter gradients.
%
% Input: 
% m
% grad_m
% usr_par : auxiliary user defined parameters (optional)
%
% Output:
% grad_parameters
%
% See also MAP_PARAMETERS_TO_M and MAP_GRADPARAMETERS_TO_GRADM.


[~, ~, nx, nz] = input_parameters();
grad_parameters = reshape( grad_m, nx, nz, 2 );


end