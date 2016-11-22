
clear all

[Lx,Lz,nx,nz] = input_parameters();
[X,Z] = define_computational_domain(Lx,Lz,nx,nz);
[width] = absorb_specs();


load ~/Desktop/hessian_H_m.mat
% load ~/Desktop/hessian_H_parameters.mat

% load ~/Desktop/paper_matlab_files/hessian/dm_structure/hessian_H_m_center_model_0.mat
% load ~/Desktop/paper_matlab_files/hessian/dm_structure/hessian_H_parameters_center_model_0.mat

% load ~/Desktop/paper_matlab_files/hessian/dm_structure/hessian_H_m_center_model_00.mat
% load ~/Desktop/paper_matlab_files/hessian/dm_structure/hessian_H_parameters_center_model_00.mat

% load ~/Desktop/paper_matlab_files/hessian/dm_structure/hessian_H_m_center_model_5.mat
% load ~/Desktop/paper_matlab_files/hessian/dm_structure/hessian_H_parameters_center_model_5.mat

% load ~/Desktop/paper_matlab_files/hessian/dm_structure/model_137/hessian_H_m_center.mat
% load ~/Desktop/paper_matlab_files/hessian/dm_structure/model_137/hessian_H_parameters_center.mat


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
set(fig1,'units','normalized','position',[.1 .3 0.5 0.4])

font_size = 8;
title_size = 12;
maker_size = 4;


subplot(1,2,1)
set( gca, 'FontSize', font_size );

hold on
% m = max(max(abs(tmp(:,:,1))));
m = 1;
mesh( X, Z, tmp(:,:,1)'/m )

offset = max(max(abs( tmp(:,:,1) )));
plot3( usr_par.network.array(:,1), usr_par.network.array(:,2), 0*usr_par.network.array(:,2) + offset, 'wv', 'MarkerSize', maker_size, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0. 0. 0.]);

plot3([width,Lx-width],[width,width],[offset,offset],'k--')
plot3([width,Lx-width],[Lz-width,Lz-width],[offset,offset],'k--')
plot3([width,width],[width,Lz-width],[offset,offset],'k--')
plot3([Lx-width,Lx-width],[width,Lz-width],[offset,offset],'k--')

colormap(cm)
c1 = 0.9 * max(max(abs(tmp(:,:,1))) );
caxis( [ -c1, c1 ] )
view([0 90])
axis square
box on

title('change of source kernel', 'FontSize', title_size)
xlabel('x [km]')
ylabel('z [km]')

xlabels = [0 1000 2000];
ylabels = [0 1000 2000];

set(gca, 'XTick', xlabels);
set(gca, 'YTick', ylabels);
% set(gca,'yaxislocation','right');

shading interp
grid on
box on
ax = gca;
ax.LineWidth = 2;
axis square
colorbar

% return




% figure
% clf

subplot(1,2,2)
set( gca, 'FontSize', font_size );

hold on
% m = max(max(abs(tmp(:,:,2))));
m = 1;
mesh( X, Z, tmp(:,:,2)'/m )

offset = 10*max(max(abs( tmp(:,:,2) )));
plot3( usr_par.network.array(:,1), usr_par.network.array(:,2), 0*usr_par.network.array(:,2) + offset, 'wv', 'MarkerSize', maker_size, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0. 0. 0.]);

plot3([width,Lx-width],[width,width],[offset,offset],'k--')
plot3([width,Lx-width],[Lz-width,Lz-width],[offset,offset],'k--')
plot3([width,width],[width,Lz-width],[offset,offset],'k--')
plot3([Lx-width,Lx-width],[width,Lz-width],[offset,offset],'k--')

colormap(cm)
c1 = 0.9 * max(max(abs(tmp(:,:,2))));
caxis( [ -c1 c1 ] )
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
colorbar


orient landscape

