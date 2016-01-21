
[Lx,Lz,nx,nz,dt,nt,~,~,~,n_basis_fct] = input_parameters();
t = -(nt-1)*dt:dt:(nt-1)*dt;


f_min = 1/7 - 0.005;
f_max = 1/7 + 0.005;


% load and filter data
data = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_0_uniform_1gaussian_homogeneous.mat');
c_data = filter_correlations( data.c_data, t, f_min, f_max );


% load inversion results
path = '~/Desktop/runs/inversion_basis_fct/source/d0_u2g_h_filtered/i10_u_h/loga_c/f2_lbfgs/';
model = dir([path 'model_*']);
model = sort_nat({model.name});


fig1 = figure;
set(fig1,'units','normalized','position',[.1 -.3 0.5 1.1])

index = 160:180;

for i = 1:length(model)
    
    clf
    plot_recordings(data.c_data(index,:),t,'vel','k',true);
    
    load([path char(model(i))])
    
    if( i==1 )
        plot_recordings(c(index,:),t,'vel','r',true);
    else
        plot_recordings(cn(index,:),t,'vel','r',true);
    end
 
    drawnow
    pause(0.5)
    
end