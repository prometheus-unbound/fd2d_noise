
function [grad_m] = map_gradparameters_to_gradm( m, grad_parameters, usr_par )
% MAP_GRADPARAMETERS_TO_GRADM function to transform the partial derivatives
% w.r.t. the physical parameters to the partial derivatives
% w.r.t. to the model variable m. This requires to apply the chain rule
% which is implicitly given by the function MAP_M_TO_PARAMETERS.
%
% Input:
% m
% grad_parameters
% usr_par : auxiliary user defined parameters (optional)
%
% Output:
% grad_m
%
% See also MAP_M_TO_PARAMETERS.


temp = 0.0 * grad_parameters;

[absbound] = init_absbound();
[ix,iz] = find( absbound == 1, 1, 'first' );

temp(ix:end-ix+1, iz:end-iz+1, 1:end-1) = imfilter( grad_parameters(ix:end-ix+1, iz:end-iz+1, 1:end-1), usr_par.kernel.imfilter.source, 'circular' );
temp(ix:end-ix+1, iz:end-iz+1, end) = imfilter( grad_parameters(ix:end-ix+1, iz:end-iz+1, end), usr_par.kernel.imfilter.structure, 'circular' );

temp(:,:,end) = usr_par.initial.mu_0 * temp(:,:,end);

grad_m = reshape(temp, [], 1);


end