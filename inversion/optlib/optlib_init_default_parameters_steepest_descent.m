
function [options] = optlib_init_default_parameters_steepest_descent(options)
    
    if (~isfield (options, 'init_step_length'))
       options.init_step_length = 1.0; 
    end

    if (~isfield (options, 'stepsize'))
       options.stepsize = 1.0; 
    end
    
    if (~isfield (options, 'verbose') )
        options.verbose = false;
    end

    if (~isfield (options, 'output_file') )
        options.output_file = 'iterations_steepest_descent.tab';
    end

    if (~isfield (options, 'max_iterations') )
        options.max_iterations = 100;
    end

    if (~isfield (options, 'tolerance') )
        options.tolerance = 1e-3;
    end

end
