
function [H_m] = map_Hparameters_to_Hm( H_parameters, usr_par )
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


temp = 0.0 * H_parameters;

[absbound] = init_absbound();
[ix,iz] = find( absbound == 1, 1, 'first' );


% temp(ix:end-ix+1, iz:end-iz+1, :) = imfilter( H_parameters(ix:end-ix+1, iz:end-iz+1, :), usr_par.kernel.imfilter.source, 'circular' );


[Lx, Lz] = input_parameters();
[~, ~, x, z] = define_computational_domain(Lx, Lz, usr_par.config.nx, usr_par.config.nz);
for i = 1:size(temp,3)-1
    temp(ix:end-ix+1, iz:end-iz+1, i) = gaussblur2d( temp(ix:end-ix+1, iz:end-iz+1, i), x(ix:end-ix+1), z(iz:end-iz+1), usr_par.kernel.sigma.source );
end
temp(ix:end-ix+1, iz:end-iz+1, end) = gaussblur2d( temp(ix:end-ix+1, iz:end-iz+1, end), x(ix:end-ix+1), z(iz:end-iz+1), usr_par.kernel.sigma.structure );
     

temp(:,:,end) = usr_par.initial.mu_0 * temp(:,:,end);

H_m = reshape(temp, [], 1);


end