
%==========================================================================
% absorbing boundaries
%
% [width, absorb_left, absorb_right, absorb_top, absorb_bottom] = absorb_specs()
%
% output:
%--------
% width: width of the boundary layer [m]
% absorb_left: absorb waves on the left boundary
% absorb_right: absorb waves on the right boundary
% absorb_top: absorb waves on the top boundary
% absorb_bottom: absorb waves on the bottom boundary
%
%==========================================================================


function [width, absorb_left, absorb_right, absorb_top, absorb_bottom] = absorb_specs()
     
    width = 50000;              % in meters

    absorb_left = 1;            % 0 or 1
    absorb_right = 1;           % 0 or 1
    absorb_top = 1;             % 0 or 1
    absorb_bottom = 1;          % 0 or 1


end