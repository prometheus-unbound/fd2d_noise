
clear all

[Lx,Lz,nx,nz] = input_parameters();
[X,Z] = define_computational_domain(Lx,Lz,nx,nz);
[width] = absorb_specs();

X = X - 1e6;
Z = Z - 1e6;



% load ~/Desktop/runs/hessian_kernels/hessian_joint_dsource_center_wd.mat
% load ~/Desktop/runs/hessian_kernels/hessian_joint_dsource_left_wd_sym.mat

load ~/Desktop/runs/hessian_kernels/hessian_joint_dstructure_center_wd_sym.mat
% load ~/Desktop/runs/hessian_kernels/hessian_joint_dstructure_left_wd_sym.mat


% load ~/Desktop/runs/hessian_kernels/hessian_joint_dsource_center_loga_sym_freqsamp_1.mat
% load ~/Desktop/runs/hessian_kernels/hessian_joint_dsource_center_loga_sym.mat



tmp = reshape( Hdm, nx, nz, 2 );
cm = cbrewer('div','RdBu',120,'PCHIP');
X = X/1000;  Z = Z/1000;  width = width/1000;
Lx = Lx/1000;  Lz = Lz/1000;
usr_par.network.array = usr_par.network.array / 1000;


figure
clf

set( gca, 'FontSize', 18 );

c1 = 0.2;

hold on
m = max(max(abs(tmp(:,:,1))));
mesh( X, Z, tmp(:,:,1)'/m )

offset = max(max(abs( tmp(:,:,1) )));
plot3( usr_par.network.array(:,1), usr_par.network.array(:,2),  0*usr_par.network.array(:,2) + offset, 'gx', 'MarkerSize', 8 )

plot3([width,Lx-width],[width,width],[offset,offset],'k--')
plot3([width,Lx-width],[Lz-width,Lz-width],[offset,offset],'k--')
plot3([width,width],[width,Lz-width],[offset,offset],'k--')
plot3([Lx-width,Lx-width],[width,Lz-width],[offset,offset],'k--')

colormap(cm)
caxis( [ -c1, c1 ] )
view([0 90])
axis image
box on

xlabel('x [km]')
ylabel('z [km]')
xlabels = [-500 0 500];
ylabels = [-200 0 200];
xlim([-749 749])
ylim([-200 200])
set(gca, 'XTick', xlabels);
set(gca, 'YTick', ylabels);
% set(gca,'yaxislocation','right');
title('change of source kernel', 'FontSize', 20)



return


set( gca, 'FontSize', 18 );
c1 = 0.2;

hold on
m = max(max(abs(tmp(:,:,2))));
mesh( X, Z, tmp(:,:,2)'/m )

offset = max(max(abs( tmp(:,:,2) )));
plot3( usr_par.network.array(:,1), usr_par.network.array(:,2),  0*usr_par.network.array(:,2) + offset, 'gx', 'MarkerSize', 8 )

plot3([width,Lx-width],[width,width],[offset,offset],'k--')
plot3([width,Lx-width],[Lz-width,Lz-width],[offset,offset],'k--')
plot3([width,width],[width,Lz-width],[offset,offset],'k--')
plot3([Lx-width,Lx-width],[width,Lz-width],[offset,offset],'k--')

colormap(cm)
caxis( [ -c1 c1 ] )
view([0 90])
axis image
box on

xlabel('x [km]')
% ylabel('z [km]')
xlabels = [-500 0 500];
ylabels = [-200 0 200];
xlim([-749 749])
ylim([-200 200])
set(gca, 'XTick', xlabels);
set(gca, 'YTick', []);
% set(gca,'yaxislocation','right');
title('change of structure kernel', 'FontSize', 20)
colorbar

