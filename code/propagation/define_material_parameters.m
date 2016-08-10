function [structure] = define_material_parameters( make_plots )

%==========================================================================
% generate material parameters mu [N/m^2] and rho [kg/m^3]
%==========================================================================
    

    % define make_plots if not specified
    if( nargin < 1 )
        make_plots = 'yes';
    end
    
    
    [~, ~, nx, nz, ~, ~, ~, model_type] = input_parameters();
    
    
    if( model_type == 1 )

        structure.rho = 3000.0 * ones(nx, nz);
        structure.mu = 4.8e10 * ones(nx, nz);    

    end

    
    if( strcmp(make_plots, 'yes') )

        array = [];
        plot_models( sqrt(structure.mu./structure.rho), [], array, [0 0 0 0] );

    end

    
end


