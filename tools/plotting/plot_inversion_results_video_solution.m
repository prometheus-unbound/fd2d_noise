
% clear all
close all
clc

path = '~/Desktop/solution.mat';
load(path);

objective = zeros(length(usr_par.model),1);
normg = zeros(length(usr_par.model),1);
max_grad = -inf;
for i = 1:length(usr_par.model)
    objective(i) = usr_par.model(i).objective;
    normg(i) = usr_par.model(i).normg;
    
    max_grad = max( max_grad, max(abs(usr_par.model(i).gradient)) );
end

fig1 = figure(1);
set(fig1,'units','normalized','position',[0.1 0.3 0.7 0.5])
subplot(2,2,3)
semilogy(1:length(usr_par.model),objective,'kx--')
title('objective function')
hold on

subplot(2,2,4)
semilogy(1:length(usr_par.model),normg,'bx--')
title('norm of gradient')
hold on


% get configuration
[~,~,nx,nz,dt,nt,~,~,~,n_basis_fct] = input_parameters();
if(n_basis_fct == 0)
    n_basis_fct = 1;
end


for i = 1:length(usr_par.model)



m_parameters = reshape( map_m_to_parameters(usr_par.model(i).m, usr_par), nx, nz, n_basis_fct );
% [dist_true,~,clim] = make_noise_source('yes');


% cm = cbrewer('div','RdBu',100,'PCHIP');
% m = max(max(max(dist_inverted)));
% clim = [-1*m 1*m];
if( ~exist('cm','var') )
    cm = [];
end

array = [];

if( exist('clim','var') )
    plot_models( m_parameters, array, cm, [clim(1) clim(2)] );
else
    plot_models( m_parameters, array, cm, [] );
end


figure(1)
hold on
if(i>1)
    subplot(2,2,3)
    semilogy(i-1,usr_par.model(i-1).objective,'kx','MarkerSize',10)
    
    subplot(2,2,4)
    semilogy(i-1,usr_par.model(i-1).normg,'bx','MarkerSize',10)
end

subplot(2,2,3)
semilogy(i,usr_par.model(i).objective,'rx','MarkerSize',10)

subplot(2,2,4)
semilogy(i,usr_par.model(i).normg,'rx','MarkerSize',10)


subplot(2,2,2)
pcolor((reshape( usr_par.model(i).gradient, nx, nz))')
shading interp
view([0 90])
% zlim([-0.5*max_grad 0.5*max_grad])
caxis([-0.005*max_grad 0.005*max_grad])

pause(0.5)

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
% plot correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_0_uniform_1gaussian_homogeneous.mat');
% initial = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_0_uniform1_homogeneous.mat');
% t = -(nt-1)*dt:dt:(nt-1)*dt;
% 
% figure
% plot_recordings(data.c_data,t,'vel','k',true);
% plot_recordings(initial.c_data,t,'vel','r',true);
% plot_recordings(model_final.cn,t,'vel','b',true);


