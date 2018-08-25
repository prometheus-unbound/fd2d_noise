

clear all
clc
close all


load random_0.07_norm.mat


[Lx, Lz, nx, nz] = input_parameters();
[X,Z,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);
[width] = absorb_specs();


fig1 = figure(1);
set(fig1,'units','normalized','position',[.1 .6 0.2 0.3])
hold on
mesh( X, Z, signal' )
level = [3 3];
plot3([width,Lx-width],[width,width],level,'k--')
plot3([width,Lx-width],[Lz-width,Lz-width],level,'k--')
plot3([width,width],[width,Lz-width],level,'k--')
plot3([Lx-width,Lx-width],[width,Lz-width],level,'k--')

caxis([-max(max(abs(signal))) max(max(abs(signal)))])
cm = cbrewer('div','RdBu',120,'PCHIP');
colormap(cm);
view([0 90])
axis square
colorbar