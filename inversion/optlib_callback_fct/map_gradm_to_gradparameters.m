
function [grad_parameters] = map_gradm_to_gradparameters(m, grad_m, usr_par)
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


[~,~,nx,nz,~,~,~,~,~,n_basis_fct] = input_parameters();

if( n_basis_fct == 0 )
    n_basis_fct = 1;
end

grad_parameters = reshape( grad_m, nx, nz, n_basis_fct + 1 );
grad_parameters(:,:,end) = grad_parameters(:,:,end) / 4.8e10;


end