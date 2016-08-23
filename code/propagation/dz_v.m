
%==========================================================================
% compute z-derivative of velocity
%
% [out] = dz_v( v, dz, nx, nz, order )
%
% input:
%-------
% v: velocity
% dz: grid spacing in z-direction
% nx, nz: number of grid points in x- and z-direction
% order: finite-difference order (2, 4, 6 or 8)
%
% output:
%--------
% velocity derivative in z-direction
%
%==========================================================================


function [out] = dz_v(v, dz, nx, nz, order)


    out = zeros(nx, nz - 1);

    if (order == 2)

        for j = 1:nz - 1
            out(:, j) = (v(:, j + 1) - v(:, j)) / dz;
        end

    elseif (order == 4)

        for j = 2:nz - 2
            out(:, j) = 9 * (v(:, j + 1) - v(:, j)) / (8 * dz) - (v(:, j + 2) - v(:, j - 1)) / (24 * dz);
        end

    end


end
