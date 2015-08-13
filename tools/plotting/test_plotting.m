

% get configuration
[f_sample,n_sample] = input_interferometry();
[Lx,Lz,nx,nz,~,~,~,model_type,source_type,n_basis_fct] = input_parameters();
[X,Z] = define_computational_domain(Lx,Lz,nx,nz);
[width] = absorb_specs();


load('~/Desktop/model_9.mat','xn')
noise_source_distribution = reshape(xn,nx,nz,n_sample);


fig1 = figure;
set(fig1,'units','normalized','position',[.1 -.3 0.5 1.1])
hold on

fudge_factor = 10;
int_limits = integration_limits(n_sample,n_basis_fct);

for ib = 8%1:n_basis_fct
    
    noise_dist_basis = noise_source_distribution(:,:,int_limits(ib,1));
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

if( exist('c_lim','var') )
    caxis(c_lim)
end
    
view([0 90])
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
