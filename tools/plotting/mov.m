
% clear all
close all


writerObj = VideoWriter('~/Desktop/test','MPEG-4');
writerObj.FrameRate = 6;
open(writerObj);

load ~/Desktop/C_out_7.mat
load ~/Desktop/runs/complex_strong/array_16_ref.mat


[Lx,Lz,nx,nz,~,~,~,model_type] = input_parameters();
[X,Z] = define_computational_domain(Lx,Lz,nx,nz);
[width] = absorb_specs();

X = X/1000;  Z = Z/1000;
width = width/1000;
array = array/1000;

fig1 = figure(1);
set(fig1,'units','normalized','position',[.1 .3 0.4 0.6])

cm = cbrewer('div','RdBu',120,'PCHIP');
offset = max(max(max(abs(x))));
c = 0.4;

for i = 50:2:size( x, 3 )-50
   
   clf
   set(gca,'FontSize',18);
   hold on
   
   mesh( X, Z, x(:,:,i)' )
   plot3( array(:,1), array(:,2), offset + 0*array(:,1), 'o' )
   plot3( array(7,1), array(7,2), offset + 0*array(7,1), 'x' )
   
   view([0 90])
   cb = colorbar;
   ylabel(cb,'[m]')
   
   colormap(cm)
   caxis([ -c*max(max(abs(x(:,:,i)))) c*max(max(abs(x(:,:,i))))])
   
   xlabel('x [km]')
   ylabel('z [km]')
   xlabels = [0 1000 2000];
   ylabels = [0 1000 2000];
   set(gca, 'XTick', xlabels);
   set(gca, 'YTick', ylabels);
   
   box on
   
   axis square
   drawnow
   
    
    M = getframe(fig1);
    writeVideo(writerObj,M); 
    
end



close(writerObj);