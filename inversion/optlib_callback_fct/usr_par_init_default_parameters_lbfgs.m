
function [usr_par] = usr_par_init_default_parameters_lbfgs(usr_par)
     

    if( ~isfield( usr_par, 'type') )
        usr_par.type = 'source';
    end
    
    
    if( ~isfield( usr_par, 'use_mex') )
        usr_par.use_mex = 'no';
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
