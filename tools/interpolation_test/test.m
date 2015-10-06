
clear all
close all
clc


load backup.mat
backup = real( backup );

cm = cbrewer('div','RdBu',100,'PCHIP');

min_backup = min(min(min( backup )));
max_backup = max(max(max( backup )));
yy = zeros(size(backup,1),size(backup,2),4);

for i = 1:size(backup,3)
    
    if( mod(i,5) ~= 1 )
        
        if( j == 1 )
            for m = 1:size(backup,1)
                for l = 1:size(backup,2)
                    yy(m,l,1:4) = spline([plot_i plot_i+5], [backup(m,l,plot_i) backup(m,l,plot_i+5)], i:(i+3) );
                end
            end
        end

        mesh( backup(:,:,i) - yy(:,:,j) )
        
%         view([0 0])
%         caxis([ -max_backup max_backup ])
        colormap(cm)
%         zlim([ -max_backup max_backup ])
        
        j = j+1;
        drawnow        
        pause(0.5)
        continue
        
    end
    
    
    mesh( 0*backup(:,:,i) )
    plot_i = i;
    j = 1;
    
%     view([0 0])
%     caxis([ -max_backup max_backup ])
    colormap(cm)
%     zlim([ -max_backup max_backup ])
    
    drawnow
    pause(0.5)
    
end

