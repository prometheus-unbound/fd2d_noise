
%==========================================================================
% user input
%==========================================================================

% small test array for gradient validation
% array = zeros(2,2);
% array(1,1) = 2.5e4;
% array(2,1) = 3.5e4;
% array(:,2) = 3.0e4;

% design array [m]
% x-components: array(:,1)
% z-components: array(:,2)
array = zeros(2, 2);
array(1, 1) = 1.4e5;
array(2, 1) = 2.6e5;
array(1:2, 2) = 2.0e5;

% array(3:4, 1) = 2.0e5;
% array(3, 2) = 1.4e5;
% array(4, 2) = 2.6e5;

% nr_x = 4;
% nr_z = 4;
% array = zeros(nr_x*nr_z,2);
% for i = 1:nr_x
%     for j = 1:nr_z        
%         array( (i-1)*nr_z + j, 1 ) = 0.8e5 + ( i-1 ) * 0.8e5;
%         array( (i-1)*nr_z + j, 2 ) = 0.8e5 + ( j-1 ) * 0.8e5;
%     end
% end

% select receivers that will be reference stations
ref_stat = array(1,:);
% ref_stat = array([6 7 10 11],:);


%==========================================================================
% calculate correlations
%==========================================================================

%- check path -------------------------------------------------------------
fd2d_path();


%- get configuration and set up time vector -------------------------------
[~, ~, nx, nz, dt, nt, ~, model_type, source_type, ~, make_plots] = input_parameters();
t = - (nt - 1) * dt:dt:(nt - 1) * dt;
nt = length(t);


%- get source and material ------------------------------------------------
noise_source = make_noise_source('no');
structure = define_material_parameters('no');


%- for gradient test ------------------------------------------------------
% noise_source.distribution = noise_source.distribution + rand(nx,nz);
% structure.mu = structure.mu + 1e9;


%- plot model with array configuration ------------------------------------
if (strcmp(make_plots, 'yes'))
    plot_models(sqrt(structure.mu ./ structure.rho), ...
        noise_source.distribution, array, [0, 0, 0, 0]);
end


%- loop over reference stations -------------------------------------------
n_ref = size(ref_stat, 1);
n_rec = size(array, 1) - 1;
correlations = zeros(n_ref, n_rec, nt);

tic
for i_ref = 1:n_ref
    
    src = ref_stat(i_ref,:);
    rec = array(~ismember(array, src, 'rows'),:);
    
    if (~exist(filename('G_fft', i_ref), 'file'))
        fprintf('ref %i: calculate Green function\n', i_ref)
        G_fft = run1_forward_green(structure, src, 0);
        parsave(filename('G_fft', i_ref), G_fft)
    else
        fprintf('ref %i: load pre-computed Green function\n', i_ref)
        G_fft = parload(filename('G_fft', i_ref));
    end
    
    fprintf('ref %i: calculate correlations\n', i_ref)
    [correlations(i_ref,:,:)] = run2_forward_correlation(structure, noise_source, G_fft, src, rec, 0);
    
    fprintf('ref %i: done\n', i_ref)
    
end
toc


%- save array and data for inversion --------------------------------------
save(filename('array', n_ref), 'array', 'ref_stat')
save(filename('correlations', n_ref), 'correlations', 't')


%- plot data --------------------------------------------------------------
% if (strcmp(make_plots, 'yes'))
    
    fig = figure;
    set(fig, 'units', 'normalized', 'position', [0.1, 0.3, 0.6, 0.5])
    plot_recordings(correlations, t, 'k', true);
    
% end


% if( strcmp(make_plots,'yes') )
% 
%     fig = figure;
%     set(fig,'units','normalized','position',[0.1 0.3 0.6 0.5])
%     list = dir([fd2d_path 'output' filesep 'correlations_nref*']);
%     handle = [];
%     legend_string = [];
%     colors = hsv(size(list,1));
%     for i = 1:size(list,1)
% 
%         tmp = load( list(i).name );
%         handle(end + 1,:) = plot_recordings(tmp.correlations, t, colors(i,:), false);
%         legend_string{end+1} = list(i).name;
% 
%     end
% 
%     legend(handle,legend_string,'Interpreter','none','Location','best')
% 
% end


