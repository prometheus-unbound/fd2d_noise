
function [Lx, Lz, nx, nz, dt, nt, order, model_type, source_type, store_fwd_nth, make_plots, plot_nth] = input_parameters()


%==========================================================================
% set basic simulation parameters
%==========================================================================

% tiny setup for gradient validation
% Lx = 6.0e4;                 % model extension in x-direction [m]
% Lz = 6.0e4;                 % model extension in y-direction [m]
% nx = 50;                    % grid points in x-direction
% nz = 50;                    % grid points in z-direction
% 
% dt = 0.09;                  % time step [s]
% nt = 50;                    % number of iterations


% setup for kernel calculation
Lx = 4.0e5;                 % model extension in x-direction [m]
Lz = 4.0e5;                 % model extension in y-direction [m]
nx = 300;                   % grid points in x-direction
nz = 300;                   % grid points in z-direction

dt = 0.09;                  % time step [s]
nt = 900;                   % number of iterations


order = 4;                  % finite-difference order (2 or 4)


%==========================================================================
% set basic parameters for kernel calculation
%==========================================================================

store_fwd_nth = 1;           % store forward wavefield every nth time step


%==========================================================================
% model type
%==========================================================================

model_type = 1;

% 1 = homogeneous 
% 2 = ?


%==========================================================================
% source type
%==========================================================================

source_type = 'homogeneous';
% source_type = 'point';
% source_type = 'gaussian';


%==========================================================================
% plotting parameters
%==========================================================================

make_plots = 'yes';         % 'yes' or 'no'
plot_nth = 100;             % plot every nth time step


end
