function plot_noise_structure_kernels(X,Z,K_rho)

%==========================================================================
%- plot noise source kernels as a function of frequency 
%==========================================================================

%- initialisations --------------------------------------------------------
cm = cbrewer('div','RdBu',100,'PCHIP');
[Lx,Lz] = input_parameters();


%- plot source and receiver positions -------------------------------------
% figure
% hold on
%     
% for k=1:length(src_x)
%     plot(src_x(k),src_z(k),'kx')
% end
% 
% for k=1:length(rec_x)
%     plot(rec_x(k),rec_z(k),'ko')
% end


%- plot noise source kernel -----------------------------------------------
m = max(max(K_rho));
kernel = pcolor(X,Z,K_rho'/m);
shading interp
alpha(kernel,0.8)


%- axis, scaling, etc. ----------------------------------------------------
colormap(cm);
set(gca,'FontSize',20);
axis equal
xlim([0 Lx])
ylim([0 Lz])

xlabel('x [m]');
ylabel('z [m]');

% m = max(max(K_rho));
% caxis([-m m]);
caxis([-1 1]);


