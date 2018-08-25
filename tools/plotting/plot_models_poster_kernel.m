% plot_noise_sources(noise_source_distribution,array,cm_psd,clim,overlay)
%
% m_parameter
% array (optional)
% clim: colormap limits (optional)
% overlay: 'yes' or 'no', only of use for n_basis_fct == 0 (optional)
%
% optional: leave empty if not wanted

function plot_models_poster_kernel( m_parameter, n_basis_fct, array, clim, overlay, video, cm_source, cm_structure )

    if( nargin < 5 || isempty(overlay) )
        overlay = 'no';
    end
    
    if( nargin < 6 || isempty(video) )
        video = 'no';
    end

    %% configuration
    [f_sample,n_sample] = input_interferometry();
    [Lx,Lz,nx,nz,~,~,~,model_type] = input_parameters();
    [X,Z] = define_computational_domain(Lx,Lz,nx,nz);
    [width] = absorb_specs();
    
    % convert to km
    X = X / 1000 - 1000;  Lx = Lx / 1000 - 1000;
    Z = Z / 1000 - 1000;  Lz = Lz / 1000 - 1000;
    width = width/1000;
    
    
    %% colormaps
    if( nargin < 7 || isempty(cm_source) )
        % cm_source = cbrewer('div','BrBG',120,'PCHIP');
        cm_source_orig = cbrewer('div','RdBu',120,'PCHIP');
        cm_source = cm_source_orig(50:120,:);
        % cm_source = cm_source_orig(47:120,:);
        % cm_source = cm_source_orig(13:120,:);
    end
    
    if( nargin < 8 || isempty(cm_structure) )
        cm_structure = cbrewer('div','RdBu',120,'PCHIP');
    end
    
    
    
    %% open figure, set size, etc.
    font_size = 18;
    title_size = 24;
    maker_size = 10;
    
    if(strcmp(video,'no') )
        fig1 = figure;
    else
        fig1 = figure(1);
    end
    
    clf
%     set(fig1,'units','normalized','position',[.1 .3 0.5 0.4])
    set(fig1,'units','normalized','position',[0.1300 0.2873 0.7201 0.6377])

    if(strcmp(video,'yes'))
        ax1 = subplot(1,2,1);
    else
        ax1 = subplot(1,1,1);
    end
    cla
    set(ax1,'FontSize',font_size);
    hold on
    
    
    % plot array if given
    if( ~isempty(array) )
        % plot3( array(:,1)/1000, array(:,2)/1000, 7+0*array(:,2), 'd', 'MarkerFaceColor', 'k', 'MarkerSize', maker_size )
        handle(1,:) = plot3( array(:,1)/1000-1000, array(:,2)/1000-1000, 7+0*array(:,2), 'wv', 'MarkerSize', maker_size, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0. 0. 0.]);
        % legend('array')
    end
    
    
    % plot for n_basis_fct == 0, i.e. only one map, or structure model
    if( n_basis_fct == 0 )
        
        if(strcmp(overlay,'yes'))

            [mu,~] = define_material_parameters(nx,nz,model_type);
            pcolor(X,Z,(mu-4.8e10)'/max(max(abs(mu-4.8e10))))

            dist = pcolor(X, Z, sum(m_parameter,3)' / max(max( sum(m_parameter,3) )) );
            alpha(dist,0.5)
            cb = colorbar;
            ylabel(cb,'normalized for overlay')
            caxis([-1.0 1.0])

        else
            
            h(1) = surf(X, Z, m_parameter(:,:,1)' )

        end
        
%         handle(2,:) = plot3([width,Lx-width],[width,width],[2,2],'k--');
%         plot3([width,Lx-width],[Lz-width,Lz-width],[2,2],'k--')
%         plot3([width,width],[width,Lz-width],[2,2],'k--')
%         plot3([Lx-width,Lx-width],[width,Lz-width],[2,2],'k--')
%         legend(handle, 'receiver', 'absorbing boundaries')
        
    % plot maps for different frequency maps
    else

        fudge_factor = 10;

        int_limits = integration_limits(n_sample,n_basis_fct);
        for ib = 1:n_basis_fct

            h(ib) = mesh(X, Z, ib*fudge_factor + m_parameter(:,:,ib)' );
            set(h(ib),'CData',m_parameter(:,:,ib)');

            text(1e3,1e3, ib*fudge_factor + 1 + m_parameter(1,1,ib) , sprintf('%5.3f - %5.3f Hz',f_sample(int_limits(ib,1)),f_sample(int_limits(ib,2))) )
            level = [ib*fudge_factor + 0.1 + m_parameter(1,1,ib), ib*fudge_factor + 0.1 + m_parameter(1,1,ib)];
            plot3([width,Lx-width],[width,width],level,'k--')
            plot3([width,Lx-width],[Lz-width,Lz-width],level,'k--')
            plot3([width,width],[width,Lz-width],level,'k--')
            plot3([Lx-width,Lx-width],[width,Lz-width],level,'k--')

        end
       
        view([22 4])
        zlim([0 fudge_factor*(n_basis_fct+1)])
        set(gca,'ZTick',[])
        zlabel('frequency bands')

    end

    
%     cb = colorbar('EastOutside');
%     cb = colorbar('SouthOutside');
%     ylabel(cb,'[kg^2 m^{-2} s^{-2}]')
%     colormap(ax1,cm_source)
%     set(cb,'YTick',[0.9,1])
            
%     if( clim(1)~=0 || clim(2)~=0 )
%         caxis([clim(1) clim(2)]);
%         set(cb,'YTick',[clim(1) clim(2)])
%     elseif( strcmp(video, 'yes') )
%         caxis([0 7])
%         set(cb,'YTick',[0 2 4 6])
%     else
%         caxis([0 7])
%         set(cb,'YTick',[0 2 4 6])
%     end

    
    xlabels = [-600 0 600];
    ylabels = [-300 0 300];    
    set(gca, 'XTick', xlabels);
    set(gca, 'YTick', ylabels);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
%     xlabel('x [km]')
%     ylabel('z [km]')
    xlim([-750 750])
    ylim([-400 400])
    
    shading interp
    grid on
    box on
    ax = gca;
    ax.LineWidth = 2;
    daspect([max(daspect)*[1 1] 1])
    
%     title('all hyperbolas', 'FontSize', title_size)
    title('window around time lag zero', 'FontSize', title_size)
%     title('surface wave on causal branch', 'FontSize', title_size)
%     title('asymmetry measurement', 'FontSize', title_size)
    
    orient landscape
%     return
    

    %% structure plot
%     if(strcmp(video,'yes'))
%         ax2 = subplot(1,2,2);
%     else
%         ax2 = subplot(1,1,1);
%     end
%     cla
%     set(ax2,'FontSize',font_size);
%     hold on
%     
%     % plot array if given
%     if( ~isempty(array) )
%         plot3( array(:,1)/1000 - 1000, array(:,2)/1000 - 1000, 2+0*array(:,2), 'wv', 'MarkerSize', maker_size, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', [0. 0. 0.]);
%     end
%     
%     
%     h(1) = surf(X, Z, m_parameter(:,:,end)');
%     
%     
%     if( clim(3)~=0 || clim(4)~=0 )
%         caxis([clim(3) clim(4)]);
% %         set(cb,'YTick',[clim(3) clim(4)])
%     elseif( strcmp(video,'yes') )
%         caxis([3.6 4.4]*1e3)
% %         set(cb,'YTick',[3.6 4.0 4.4]*1e3,'FontSize',28)
%     else
% %         caxis([3.87 4.13]*1e3)
% %         set(cb,'YTick',[3.9 4.0 4.1]*1e3)
%     end
%     colormap(ax2, cm_structure);
% 
%     
%     xlabels = [-500 -250 0 250];
%     ylabels = [-100 0 100];    
% %     set(gca, 'XTick', xlabels);
% %     set(gca, 'YTick', ylabels);
%     set(gca, 'XTick', []);
%     set(gca, 'YTick', []);
% %     xlabel('x [km]')
% %     ylabel('z [km]')
%     xlim([-600 400])
%     ylim([-200 200])
%                 
%     shading interp   
%     
%     grid on
%     box on
%     ax = gca;
%     ax.LineWidth = 2;
%     daspect([max(daspect)*[1 1] 1])
    
    
%     title('homogeneous distribution', 'FontSize', title_size-2)
%     title('Gaussian anomaly', 'FontSize', title_size-2)
%     title('point source - in receiver line', 'FontSize', title_size-2)
%     title('point source - off receiver line', 'FontSize', title_size-2)


    
%     % gaussian
%     center_x = -0.4e3;
%     center_z = -0.1e3;
%     width = 0.75e2;
    
%     % point
%     center_x = 0.0;
%     center_z = -0.1e3;
%     width = 0.1e2;
    
%     % point_tsm
    center_x = -0.4e3;
    center_z = 0.0;
    width = 0.1e2;
    
    X_interp = interp2(X,2);
    Z_interp = interp2(Z,2);
    pattern = 0.0 * double( (X_interp-center_x).^2 + (Z_interp-center_z).^2 <= width^2 );
    pattern(pattern==0) = NaN;
    
    h(2) = surf(X_interp, Z_interp, pattern);
    
    set(h(2), 'FaceColor', 'interp')
%     alpha(h(2), 0.35)   
%     alpha(h(2), 1.0)
    alpha(h(2), 0.0)

    c_data_1 = m_parameter(:,:,1)';
    c_data_1( c_data_1 <= clim(1) ) = 0.99*clim(1);
    c_data_1( c_data_1 >= clim(2) ) = 0.99*clim(2);    
    caxis([clim(1) clim(2) + clim(2)-clim(1)])
    
%     c_data_1 = m_parameter(:,:,end)';
%     c_data_1( c_data_1 <= clim(3) ) = 0.99*clim(3);
%     c_data_1( c_data_1 >= clim(4) ) = 0.99*clim(4);
%     caxis([clim(3) clim(4) + clim(4)-clim(3)])
    
    set( h(1), 'CData', c_data_1 );
    set( h(2), 'CData', pattern );
    
    cm = [cm_structure; cbrewer('seq','YlGn',120,'PCHIP')];
    colormap(ax1, cm);
%     colormap(ax2, cm);
    cb = colorbar;
    set(cb,'Limits',[clim(1) clim(2)])
    set(cb,'YTick',[clim(1) clim(2)])
%     set(cb,'Limits',[clim(3) clim(4)])
%     set(cb,'YTick',[clim(3) clim(4)])
    
    set(cb, 'visible', 'on')
    
    shading interp
    orient landscape
    
end

