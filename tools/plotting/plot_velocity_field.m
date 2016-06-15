
if (mod(n,plot_nt)==0)
    
    hold on
    
    %- plot velocity field ------------------------------------------------ 
    pcolor(X,Z,v');
    
    [width,absorb_left,absorb_right,absorb_top,absorb_bottom] = absorb_specs();
    plot([width,Lx-width],[width,width],'k--')
    plot([width,Lx-width],[Lz-width,Lz-width],'k--')
    plot([width,width],[width,Lz-width],'k--')
    plot([Lx-width,Lx-width],[width,Lz-width],'k--')
    
    set(gca,'FontSize',20);
    axis image
    
    
    %- scale, label, etc ... ----------------------------------------------   
    if (n<0.8*length(t))
        scale=max(max(abs(v)));
    end
   
    caxis([-scale scale]);
    colorbar
    colormap(cm);
    shading interp
    xlabel('x [m]');
    ylabel('z [m]');
    
    
    %- record movie -------------------------------------------------------
    % if exist('movie_index','var')
    %     movie_index=movie_index+1;
    % else
    %     movie_index=1;
    % end
    % 
    % M(movie_index) = getframe(gcf);
        
    
    %- finish -------------------------------------------------------------
    drawnow
    hold off
    clf
            
end
