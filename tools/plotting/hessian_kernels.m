
clear all

[Lx,Lz,nx,nz] = input_parameters();
[X,Z] = define_computational_domain(Lx,Lz,nx,nz);
[width] = absorb_specs();


% load ~/Desktop/hessian_H_m.mat
% load ~/Desktop/hessian_H_parameters_loga.mat
% load ~/Desktop/hessian_H_parameters_ampdiff.mat
% load ~/Desktop/hessian_H_parameters_ampdiff_wd.mat

% load ~/Diss/Paper/phd_paper_1/matlab_files/hessian/dm_structure/hessian_H_m_center_model_0.mat
% load ~/Diss/Paper/phd_paper_1/matlab_files/hessian/dm_structure/hessian_H_parameters_center_model_0.mat

% load ~/Diss/Paper/phd_paper_1/matlab_files/hessian/dm_structure/hessian_H_m_center_model_00.mat
% load ~/Diss/Paper/phd_paper_1/matlab_files/hessian/dm_structure/hessian_H_parameters_center_model_00.mat

% load ~/Diss/Paper/phd_paper_1/matlab_files/hessian/dm_structure/hessian_H_m_center_model_5.mat
load ~/Diss/Paper/phd_paper_1/matlab_files/hessian/dm_structure/hessian_H_parameters_center_model_5.mat

% load ~/Diss/Paper/phd_paper_1/matlab_files/hessian/dm_structure/model_137/hessian_H_m_center.mat
% load ~/Diss/Paper/phd_paper_1/matlab_files/hessian/dm_structure/model_137/hessian_H_parameters_center.mat


% load ~/Desktop/paper_matlab_files/hessian/dm_source/hessian_H_m_lefts_model_137.mat
% load ~/Desktop/paper_matlab_files/hessian/dm_source/hessian_H_parameters_lefts_model_137.mat


% usr_par.kernel.sigma.source = [1e-3 1e-3];
% usr_par.kernel.sigma.structure = usr_par.kernel.sigma.source;

if( exist('H_m','var') )
    tmp = reshape( H_m, nx, nz, 2 );
elseif( exist('H_parameters','var') )
    tmp = H_parameters;
else
    tmp = map_gradm_to_gradparameters( 0, Hdm, usr_par );
end

cm = cbrewer('div','RdBu',120,'PCHIP');
X = X/1000;  Z = Z/1000;  width = width/1000;
Lx = Lx/1000;  Lz = Lz/1000;
usr_par.network.array = usr_par.network.array / 1000;


fig1 = figure;
% set(fig1,'units','normalized','position',[.1 .3 0.5 0.4])
set(fig1,'units','normalized','position',[0.1300 0.2873 0.7201 0.6377])

font_size = 18;
title_size = 24;
maker_size = 8;


% subplot(1,1,1)
% % subplot(1,2,1)
% set( gca, 'FontSize', font_size );
% 
% hold on
% % m = max(max(abs(tmp(:,:,1))));
% m = 1;
% mesh( X, Z, tmp(:,:,1)'/m )
% 
% 
% % offset = max(max(abs( tmp(:,:,1) )));
% offset = 10.0;
% 
% dm_parameters_2 = inf * tmp;
% dm_parameters_2( 382:386, 291:295, 1 ) = offset;
% test = mesh( X, Z, dm_parameters_2(:,:,1)', 10000*dm_parameters_2(:,:,1)');
% 
% plot3( usr_par.network.array(:,1), usr_par.network.array(:,2), 0*usr_par.network.array(:,2) + offset, 'wv', 'MarkerSize', maker_size, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0. 0. 0.]);
% 
% plot3([width,Lx-width],[width,width],[offset,offset],'k--')
% plot3([width,Lx-width],[Lz-width,Lz-width],[offset,offset],'k--')
% plot3([width,width],[width,Lz-width],[offset,offset],'k--')
% plot3([Lx-width,Lx-width],[width,Lz-width],[offset,offset],'k--')
% 
% colormap(cm)
% c1 = 1.0 * max(max(abs(tmp(:,:,1))) );
% % c1 = 1.0;
% caxis( [ -c1, c1 ] )
% % caxis( [ -1, 1 ]*1e-6 )
% % caxis( [ -5, 5 ]*1e-5 )
% view([0 90])
% axis square
% box on
% 
% title('change of source kernel', 'FontSize', title_size)
% xlabel('x [km]')
% ylabel('z [km]')
% 
% xlabels = [0 1000 2000];
% ylabels = [0 1000 2000];
% 
% set(gca, 'XTick', xlabels);
% set(gca, 'YTick', ylabels);
% % set(gca,'yaxislocation','right');
% 
% shading interp
% grid on
% box on
% ax = gca;
% ax.LineWidth = 2;
% axis square
% cb = colorbar;
% % set(cb,'YTick',[-1 0 1]*1e-6)
% set(cb,'YTick',[-5 0 5]*1e-5)
% % set(cb,'YTick',[-0.8 0 0.8])
% % set(cb,'visible', 'off')
% 
% set(test,'Edgecolor',[0 0.8 0])
% orient landscape
% return




% figure
% clf
subplot(1,1,1)
% subplot(1,2,2)
set( gca, 'FontSize', font_size );

hold on
% m = max(max(abs(tmp(:,:,2))));
m = 1;
mesh( X, Z, tmp(:,:,2)'/m )

% offset = 10*max(max(abs( tmp(:,:,2) )));
offset = 10.0;


dm_parameters = inf * tmp;
dm_parameters( 382:386, 291:295, 2 ) = offset;
test = mesh( X, Z, dm_parameters(:,:,2)', 10000*dm_parameters(:,:,2)');



plot3( usr_par.network.array(:,1), usr_par.network.array(:,2), 0*usr_par.network.array(:,2) + offset, 'wv', 'MarkerSize', maker_size, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0. 0. 0.]);

plot3([width,Lx-width],[width,width],[offset,offset],'k--')
plot3([width,Lx-width],[Lz-width,Lz-width],[offset,offset],'k--')
plot3([width,width],[width,Lz-width],[offset,offset],'k--')
plot3([Lx-width,Lx-width],[width,Lz-width],[offset,offset],'k--')

colormap(cm)
c1 = 1.0 * max(max(abs(tmp(:,:,2))));
% c1 = 1.0;
caxis( [ -c1 c1 ] )
% caxis( [ -1 1 ]*1e-15 )
% caxis( [ -2 2 ]*1e-13 )
view([0 90])
axis square
box on

title('change of structure kernel', 'FontSize', title_size)
xlabel('x [km]')
% ylabel('z [km]')

xlabels = [0 1000 2000];
ylabels = [0 1000 2000];
xlim([0 2000])
ylim([0 2000])

set(gca, 'XTick', xlabels);
set(gca, 'YTick', []);
% set(gca,'yaxislocation','right');

shading interp
grid on
box on
ax = gca;
ax.LineWidth = 2;
axis square
cb = colorbar;
% set(cb,'YTick',[-1 0 1]*1e-3)
% set(cb,'YTick',[-2 0 2]*1e-13)
set(cb,'YTick',[-1 0 1]*1e-15)
% set(cb,'YTick',[-0.8 0 0.8])
set(test,'Edgecolor',[0 0.9 0])


orient landscape

