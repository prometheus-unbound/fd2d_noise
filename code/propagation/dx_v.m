
%==========================================================================
% compute x-derivative of velocity
%
% [out] = dx_v( v, dx, nx, nz, order )
%
% input:
%-------
% v: velocity
% dx: grid spacing in x-direction
% nx, nz: number of grid points in x- and z-direction
% order: finite-difference order (2, 4, 6 or 8)
%
% output:
%--------
% velocity derivative in x-direction
%
%==========================================================================


function [out] = dx_v(v, dx, nx, nz, order)


    out = zeros(nx - 1, nz);

    if (order == 2)

        for i = 1:nx - 1
            out(i,:) = (v(i + 1,:) - v(i,:)) / dx;
        end

    elseif (order == 4)

        for i = 2:nx - 2
            out(i,:) = 9 * (v(i + 1,:) - v(i,:)) / (8 * dx) - (v(i + 2,:) - v(i - 1,:)) / (24 * dx);
        end

    end


end
