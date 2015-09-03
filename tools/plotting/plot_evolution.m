
[Lx,Lz,nx,nz,dt,nt,~,~,~,n_basis_fct] = input_parameters();

t = -(nt-1)*dt:dt:(nt-1)*dt;
data = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_0_uniform_1gaussian_homogeneous.mat');

path = '~/Desktop/correlations/';
correlations = dir([path 'correlations_*']);
correlations = {correlations.name};

% path = '~/Desktop/runs/inversion_basis_fct/source/d0_u1g_h/i10_u_h/loga_c/';
% correlations = dir([path 'model_*']);
% correlations = sort_nat({correlations.name});

fig1 = figure;
set(fig1,'units','normalized','position',[.1 -.3 0.5 1.1])

index = 160:160;

for i = 1:length(correlations)
    
    clf
    plot_recordings(data.c_data(index,:),t,'vel','k',true);
    
    load([path char(correlations(i))],'c_all')
    plot_recordings(c_all(index,:),t,'vel','r',true);
 
    drawnow
    pause(0.1)
    
end