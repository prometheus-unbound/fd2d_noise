
function [usr_par] = usr_par_init_default_parameters_lbfgs(usr_par)
     

    if( ~isfield( usr_par, 'type') )
        usr_par.type = 'source';
    end
    
    
    if( ~isfield( usr_par, 'use_mex') )
        usr_par.use_mex = 'no';
    end
    
    
    if( strcmp( usr_par.use_mex, 'no' ) )
        ! rm ../code_mex_functions/run*
        ! cp ../code/run1_forward_green.m ../code/mex_functions/run1_forward_green_mex.m
        ! cp ../code/run2_forward_correlation.m ../code/mex_functions/run2_forward_correlation_mex.m
        ! cp ../code/run3_adjoint.m ../code/mex_functions/run3_adjoint_mex.m
    end
    
    
    if( ~isfield( usr_par, 'measurement') )
        usr_par.measurement = 'waveform_difference';
    end

    
    if( isfield(usr_par,'kernel') )
        if( ~isfield( usr_par.kernel, 'imfilter') ) 
            usr_par.kernel.imfilter = fspecial('gaussian',[1 1], 1);
        end
    else
        usr_par.kernel.imfilter = fspecial('gaussian',[1 1], 1);
    end
    
    
    if( ~isfield( usr_par, 'network') )
        error('\nload array!\n')
    end
    
    
    if( ~isfield( usr_par, 'data') )
        error('\nload data!\n')
    end
    

end
