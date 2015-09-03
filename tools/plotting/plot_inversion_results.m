
% c_lim = get(gca,'CLim');

% get configuration
[f_sample,n_sample] = input_interferometry();
[Lx,Lz,nx,nz,dt,nt,~,~,~,n_basis_fct] = input_parameters();
[X,Z] = define_computational_domain(Lx,Lz,nx,nz);
[width] = absorb_specs();

% path = '~/Desktop/runs/inversion_basis_fct/source/d0_u2g_h_filtered/i10_u_h/loga_c/f2_lbfgs_unfiltered/';
path = '~/Desktop/runs/inversion_basis_fct/source/d0_u2g_h_filtered/i10_u_h/loga_c/f2_lbfgs/';
% path = '~/Desktop/runs/inversion_basis_fct/source/d0_u2g_h/i10_u_h/loga_c/';

n_models = length( dir([path 'model_*']) );
model_final = load([path 'model_' num2str(n_models-1) '.mat']);
noise_source_distribution = reshape(model_final.xn,nx,nz,n_basis_fct);


fig1 = figure;
set(fig1,'units','normalized','position',[.1 -.3 0.5 1.1])
hold on

fudge_factor = 10;
int_limits = integration_limits(n_sample,n_basis_fct);

for ib = 1:n_basis_fct
    
    noise_dist_basis = noise_source_distribution(:,:,ib);
    h(ib) = mesh(X, Z, ib*fudge_factor + noise_dist_basis' );
    set(h(ib),'CData',noise_dist_basis');
    
    text(1e3,1e3, ib*fudge_factor + 1 + noise_dist_basis(1,1) , sprintf('%5.2f - %5.2f Hz',f_sample(int_limits(ib,1)),f_sample(int_limits(ib,2))) )
    level = [ib*fudge_factor + 0.1 + noise_dist_basis(1,1), ib*fudge_factor + 0.1 + noise_dist_basis(1,1)];
    plot3([width,Lx-width],[width,width],level,'k--')
    plot3([width,Lx-width],[Lz-width,Lz-width],level,'k--')
    plot3([width,width],[width,Lz-width],level,'k--')
    plot3([Lx-width,Lx-width],[width,Lz-width],level,'k--')
    
end

load cm_psd
colormap(cm_psd)
colorbar

% if( exist('c_lim','var') )
%     caxis([c_lim(1) 1/2*c_lim(2)])
% end
    
view([22 4])
% view([0 90])
zlim([0 fudge_factor*(n_basis_fct+1)])
set(gca,'ZTick',[])
zlabel('frequency bands')




load('~/Desktop/runs/inversion/data/array_16_ref.mat')
plot(array(:,1),array(:,2),'ko')

shading interp
grid on
box on
axis square
title('power-spectral density distribution');
xlabel('x [m]')
ylabel('z [m]')
xlim([0 Lx])
ylim([0 Lz])



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


