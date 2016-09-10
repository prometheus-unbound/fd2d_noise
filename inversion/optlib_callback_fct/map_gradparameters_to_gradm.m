
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


[absbound] = init_absbound();
[ix,iz] = find( absbound == 1, 1, 'first' );

[Lx, Lz, nx, nz] = input_parameters();
[~, ~, x, z] = define_computational_domain(Lx, Lz, nx, nz);

grad_parameters(ix:end-ix+1, iz:end-iz+1, 1) = gaussblur2d(grad_parameters(ix:end-ix+1, iz:end-iz+1, 1), x(ix:end-ix+1), z(iz:end-iz+1), usr_par.smoothing.sigma);
grad_parameters(ix:end-ix+1, iz:end-iz+1, 2) = gaussblur2d(grad_parameters(ix:end-ix+1, iz:end-iz+1, 2), x(ix:end-ix+1), z(iz:end-iz+1), usr_par.smoothing.sigma);

for i = 1:2
    grad_parameters(:,:,i) = double(absbound == 1) .* grad_parameters(:,:,i);
end

grad_m = reshape(grad_parameters, [], 1);


end