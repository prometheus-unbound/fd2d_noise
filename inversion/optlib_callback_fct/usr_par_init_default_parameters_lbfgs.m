
function [usr_par] = usr_par_init_default_parameters_lbfgs(usr_par)
     

    if( ~isfield( usr_par, 'type') )
        usr_par.type = 'source';
    end
    
    
    if( ~isfield( usr_par, 'use_mex') )
        usr_par.use_mex = 'no';
    end
        
    
    if( isfield( usr_par, 'measurement') )
        if( ~isfield( usr_par.measurement, 'type' ) )
            usr_par.measurement.type = 'waveform_difference';
        end
    else
        usr_par.measurement.type = 'waveform_difference';
    end
    
    
    if( isfield( usr_par, 'measurement') )
        if( ~isfield( usr_par.measurement, 'mode' ) )
            usr_par.measurement.mode = 'manual';
        end
    else
        usr_par.measurement.mode = 'manual';
    end

    
    if( isfield(usr_par,'kernel') )
        if( ~isfield( usr_par.kernel, 'imfilter') ) 
            usr_par.kernel.imfilter = fspecial('gaussian',[1 1], 1);
        end
    else
        usr_par.kernel.imfilter = fspecial('gaussian',[1 1], 1);
    end
    
    
    if( ~isfield( usr_par, 'data_independent') )
        usr_par.data_independent = 'no';
    end
    
    
    if( ~isfield( usr_par, 'network') )
        error('\nload array!\n')
    end
    
    
    if( ~isfield( usr_par, 'data') )
        if( strcmp( usr_par.data_independent, 'no' ) )
            error('\nload data!\n')
        end
    end
    

end
