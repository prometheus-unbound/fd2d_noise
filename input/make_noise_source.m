
function [ noise_source ] = make_noise_source( make_plots )


    % define make_plots if not specified
    if( nargin < 1 )
        make_plots = 'yes';
    end

    
    % get configuration
    f_sample = input_interferometry();
    [Lx, Lz, nx, nz, ~, ~, ~, ~, source_type] = input_parameters();
    [X, Z] = define_computational_domain(Lx,Lz,nx,nz);  
    
    
    %======================================================================
    % user input
    %======================================================================
    
    % specify spectrum
    f_peak = 1/9;
    bandwidth = 0.04;
    strength = 1;
    
    % point source
    if( strcmp(source_type, 'point') )
        
        x_source = 2.0e5;
        z_source = 1.7e5;
        source_width = 5e3;
        magnitude = 3.0;
        
    % gaussian anomaly
    elseif( strcmp(source_type, 'gaussian') )
        
        x_source = 1.0e5;
        z_source = 1.6e5;
        source_width = 3e4;
        magnitude = 15.0;
        
    elseif( ~strcmp(source_type, 'homogeneous')  )
       
        error('\nwrong source type!\n')
        
    end
    
   
    %======================================================================
    % set up source spectrum
    %======================================================================
    
    noise_source.spectrum = strength * exp( -( abs(f_sample)-f_peak ).^2 / bandwidth^2 );       
   
     
    %======================================================================
    % set up geographic power-spectral density distribution
    %======================================================================
    
    % homogeneous
    if( strcmp(source_type, 'homogeneous') )
        
        noise_source.distribution = ones(nx, nz);
    
    % point source
    elseif( strcmp(source_type, 'point') )
        
        noise_source.distribution = magnitude * exp( -( (X-x_source).^2 + (Z-z_source).^2 ) / source_width^2 )';
        
    % gaussian anomaly
    elseif( strcmp(source_type, 'gaussian') )
        
        noise_source.distribution = ones(nx, nz);
        noise_source.distribution = noise_source.distribution + magnitude * exp( -( (X-x_source).^2 + (Z-z_source).^2 ) / source_width^2 )';
 
    end   
    
    
    %======================================================================
    % plot noise source configuration
    %======================================================================
    
    if ( strcmp(make_plots,'yes') )
        
        figure
        set(gca,'FontSize',12)
        hold on
        grid on
        plot( f_sample, noise_source.spectrum, 'r')
        xlabel('frequency [Hz]');
        xlim([f_sample(1) f_sample(end)])
        title('spectrum for noise source','FontSize',18)
        
        if( exist([ fd2d_path() 'output' filesep 'array_1_ref.mat' ], 'file' ) )
            load([ fd2d_path() 'output' filesep 'array_1_ref.mat' ]);
        else
            array = [];
        end
        
        plot_models( [], noise_source.distribution, array, [0 0 0 90]);
        
    end
    
    
end

