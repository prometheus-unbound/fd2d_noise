
[Lx, Lz, nx, nz, dt, nt, order, model_type, source_type, n_basis_fct, fw_nth] = input_parameters();
t = -(nt-1)*dt:dt:(nt-1)*dt;


f_min = 1/7 - 0.005;
f_max = 1/7 + 0.005;


% load and filter data
% data = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_0_uniform_1gaussian_homogeneous.mat');
% c_data = filter_correlations( data.c_data, t, f_min, f_max );


% load inversion results
% path = '~/Desktop/models/';
path = '~/Desktop/models_first_run/';
models_dir = dir([path 'model_*']);
models_dir = sort_nat({models_dir.name});


fig1 = figure(1);
set(fig1,'units','normalized','position',[.1 .3 0.5 0.4])


% index = 160:180;
% for i = 1:length(model)
%     clf
%     plot_recordings(data.c_data(index,:),t,'vel','k',true);
%     
%     load([path char(model(i))])
%     
%     if( i==1 )
%         plot_recordings(c(index,:),t,'vel','r',true);
%     else
%         plot_recordings(cn(index,:),t,'vel','r',true);
%     end
%  
%     drawnow
%     pause(0.5)
% end

writerObj = VideoWriter('~/Desktop/test','MPEG-4');
writerObj.FrameRate = 6;
open(writerObj);

% load clim.mat

usr_par.network = []; usr_par.data = [];
usr_par.ring.switch = 'no';
[usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);
array = [];

for i = 1:length(models_dir)

    load([path char(models_dir(i))]);    
    
    % get true source and material
    if( usr_par.config.n_basis_fct == 0 )
        m_true = zeros( usr_par.config.nx, usr_par.config.nz, 2 );
    else
        m_true = zeros( usr_par.config.nx, usr_par.config.nz, usr_par.config.n_basis_fct+1 );
    end
    
    m_parameters = map_m_to_parameters(model.m, usr_par);
%     plot_models( m_parameters, 0, array, [0 0 clim(1) clim(2)], 'no', 'yes');
%     plot_models( m_parameters, 0, array, [0 0 4.6e10 5.0e10], 'no', 'yes');
    plot_models( m_parameters, 0, array, [0 0 4.6e10 5.0e10], 'no', 'yes');
    
    fig1 = figure(1);
    subplot(2,2,3)
    if(i==1)
        normg0 = norm(model.gradient);
    end
    semilogy(i,norm(model.gradient)/normg0,'x')
    % text(i,norm(model.gradient)/normg0,step(i))
    xlim([0 length(models_dir)])
    ylim([1e-3 1])
    % axis square
    hold on
    
    subplot(2,2,4)
    semilogy(i,misfit(i),'x')
    xlim([0 length(models_dir)])
    ylim([5 70])
    % axis square
    hold on
    
    
    fig2 = figure(2);
    if(i==1)
        normg0 = norm(model.gradient);
    end
    semilogy(i,norm(model.gradient)/normg0,'x')
    text(i,norm(model.gradient)/normg0,step(i),'FontSize',7)
    xlim([0 length(models_dir)])
    ylim([1e-3 1])
    hold on
    
    
    drawnow
    % pause(0.1)
    
    M = getframe(fig1);
    writeVideo(writerObj,M);
    
end


close(writerObj);