% plot_noise_sources(noise_source_distribution,array,cm_psd,clim,overlay)
%
% m_parameter
% array (optional)
% clim: colormap limits (optional)
% overlay: 'yes' or 'no', only of use for n_basis_fct == 0 (optional)
%
% optional: leave empty if not wanted

function plot_models( m_parameter, n_basis_fct, array, clim, overlay, video, cm_source, cm_structure )

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
    X = X / 1000;  Lx = Lx / 1000;
    Z = Z / 1000;  Lz = Lz / 1000;
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
    font_size = 8;
    title_size = 12;
    maker_size = 4;
    
    if(strcmp(video,'no') )
        fig1 = figure;
    else
        fig1 = figure(1);
    end
    
    % clf
    set(fig1,'units','normalized','position',[.1 .3 0.5 0.4])
    if(strcmp(video,'yes'))
        ax1 = subplot(1,2,1);
    else
        ax1 = subplot(1,2,1);
    end
    cla
    set(ax1,'FontSize',font_size);
    hold on
    
    
    %% plot array if given
    if( ~isempty(array) )
        % plot3( array(:,1)/1000, array(:,2)/1000, 7+0*array(:,2), 'd', 'MarkerFaceColor', 'k', 'MarkerSize', maker_size )
        plot3( array(:,1)/1000, array(:,2)/1000, 7+0*array(:,2), 'wv', 'MarkerSize', maker_size, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0. 0. 0.]);
        % legend('array')
    end
    
    
    %% plot for n_basis_fct == 0, i.e. only one map, or structure model
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
            
            mesh(X, Z, m_parameter(:,:,1)' )

        end
        
        plot3([width,Lx-width],[width,width],[2,2],'k--')
        plot3([width,Lx-width],[Lz-width,Lz-width],[2,2],'k--')
        plot3([width,width],[width,Lz-width],[2,2],'k--')
        plot3([Lx-width,Lx-width],[width,Lz-width],[2,2],'k--')

        
    %% plot maps for different frequency maps
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

    
    cb = colorbar;
    ylabel(cb,'[kg^2 m^{-2} s^{-2}]')
    colormap(ax1,cm_source)
    % set(cb,'YTick',[0.9,1])
            
    if( clim(1)~=0 || clim(2)~=0 )
        caxis([clim(1) clim(2)]);
        set(cb,'YTick',[clim(1) clim(2)])
    elseif( strcmp(video, 'yes') )
        caxis([0 7])
        set(cb,'YTick',[0 2 4 6])
    else
        % caxis([0 7])
        % set(cb,'YTick',[0 2 4 6])
    end
    
    xlabels = [0 1000 2000];
    ylabels = [0 1000 2000];
    set(gca, 'XTick', xlabels);
    set(gca, 'YTick', ylabels);
    xlabel('x [km]')
    ylabel('z [km]')
    xlim([0 Lx])
    ylim([0 Lz])
    
    shading interp
    grid on
    box on
    ax = gca;
    ax.LineWidth = 2;
    axis square
    title('source distribution', 'FontSize', title_size)
    
    
    % return
    

    %% structure plot
    if(strcmp(video,'yes'))
        ax2 = subplot(1,2,2);
    else
        ax2 = subplot(1,2,2);
    end
    cla
    set(ax2,'FontSize',font_size);
    hold on
    
    % plot array if given
    if( ~isempty(array) )
        % plot3( array(:,1)/1000, array(:,2)/1000, 6e3+0*array(:,2), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', maker_size )
        plot3( array(:,1)/1000, array(:,2)/1000, 6e3+0*array(:,2), 'wv', 'MarkerSize', maker_size, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0. 0. 0.]);
        % legend('array') 
    end
    
    
    load ../output/interferometry/array_16_ref.mat
    min_x = min(array(:,1))/1e3;
    min_z = min(array(:,2))/1e3;
    max_x = max(array(:,1))/1e3;
    max_z = max(array(:,2))/1e3;
    
    h(1) = surf(X, Z, m_parameter(:,:,end)');
    
    cb = colorbar;
    ylabel(cb,'[m/s]')
    
    
    buffer = 1.0e2;
    pattern = double( X > (min_x-buffer) & X < (max_x+buffer) ) .* double( Z > (min_z-buffer) & Z < (max_z+buffer) );
    pattern( pattern == 1 ) = NaN;
    pattern( pattern == 0 ) = 4600;    
    h(2) = mesh(X, Z, pattern );
    
    set(h(2), 'FaceColor', 'interp')
    alpha(h(2), 0.6)   

    c_data_1 = m_parameter(:,:,end)';
    c_data_1( c_data_1 > 4300 ) = 4299;
    c_data_1( c_data_1 < 3700 ) = 3700;    
    
    set( h(1), 'CData', c_data_1 );
    set( h(2), 'CData', pattern );
    
    caxis([3700 4900])
    set(cb,'Limits',[3700 4300])
    set(cb,'YTick',[3700 4000 4300])
    
    cm = [cm_structure; cbrewer('seq','Greys',120,'PCHIP')];
    colormap(ax2, cm);
    
    
%     if( clim(3)~=0 || clim(4)~=0 )
%         caxis([clim(3) clim(4)]);
%         set(cb,'YTick',[clim(3) clim(4)])
%     elseif( strcmp(video,'yes') )
%         caxis([3.6 4.4]*1e3)
%         set(cb,'YTick',[3.6 4.0 4.4]*1e3,'FontSize',28)
%     else
% %         caxis([3.87 4.13]*1e3)
% %         set(cb,'YTick',[3.9 4.0 4.1]*1e3)
%     end
%     colormap(ax2, cm_structure);

    
    xlabels = [0 1000 2000];
    ylabels = [0 1000 2000];
    set(gca, 'XTick', xlabels);
    set(gca, 'YTick', ylabels);
%     set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    xlabel('x [km]')
%     ylabel('z [km]')
    xlim([0 Lx])
    ylim([0 Lz])
        
    
    plot3([width,Lx-width],[width,width],[4.5e3,4.5e3],'k--')
    plot3([width,Lx-width],[Lz-width,Lz-width],[4.5e3,4.5e3],'k--')
    plot3([width,width],[width,Lz-width],[4.5e3,4.5e3],'k--')
    plot3([Lx-width,Lx-width],[width,Lz-width],[4.5e3,4.5e3],'k--')
    
   
%     plot3([min_x-buffer,max_x+buffer],[min_z-buffer,min_z-buffer],[4.5e3,4.5e3],'k--','LineWidth', 3)
%     plot3([min_x-buffer,max_x+buffer],[max_z+buffer,max_z+buffer],[4.5e3,4.5e3],'k--','LineWidth', 3)
%     plot3([min_x-buffer,min_x-buffer],[min_z-buffer,max_z+buffer],[4.5e3,4.5e3],'k--','LineWidth', 3)
%     plot3([max_x+buffer,max_x+buffer],[min_z-buffer,max_z+buffer],[4.5e3,4.5e3],'k--','LineWidth', 3)    

    
    shading interp   
    grid on
    box on
    ax = gca;
    ax.LineWidth = 2;
    axis square
    title('velocity model', 'FontSize', title_size)
    
        
    orient landscape
    
end

