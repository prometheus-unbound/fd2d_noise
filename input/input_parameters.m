
function [Lx,Lz,nx,nz,dt,nt,order,model_type,source_type,n_basis_fct] = input_parameters()

%==========================================================================
% set basic simulation parameters
%==========================================================================

% Lx = 4.0e5;         % model extension in x-direction [m]
% Lz = 4.0e5;         % model extension in y-direction [m]
% 
% nx = 300;           % grid points in x-direction
% nz = 300;           % grid points in z-direction
% 
% dt = 0.09;          % time step [s]
% nt = 900;           % number of iterations


Lx = 2.0e6;         % model extension in x-direction [m]
Lz = 2.0e6;         % model extension in y-direction [m]

nx = 600;           % grid points in x-direction
nz = 600;           % grid points in z-direction

dt = 0.23;          % time step [s]
nt = 1600;          % number of iterations


order=4;            % finite-difference order (2 or 4)


%==========================================================================
% model type
%==========================================================================

model_type = 999;
% model_type = 999;

% model_type = 100;

% 1=homogeneous 
% 2=homogeneous with localised density perturbation
% 3=layered medium
% 4=layered with localised density perturbation
% 5=different layered medium
% 6=vertical gradient medium in mu
% 7=another layered medium
% 666 = put source picture in ../models-folder
% 999 = IUGG structure, two pills


%==========================================================================
% source type
%==========================================================================

% source_type = 'homogeneous';
source_type = 'gaussian';

% number of frequency bands
n_basis_fct = 10;

