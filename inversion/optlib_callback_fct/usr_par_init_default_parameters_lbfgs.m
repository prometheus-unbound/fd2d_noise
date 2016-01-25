
function [usr_par] = usr_par_init_default_parameters_lbfgs(usr_par)
    

    if( ~isfield( usr_par, 'cluster') )
        usr_par.cluster = 'local';
    end


    if( ~isfield( usr_par, 'type') )
        usr_par.type = 'structure';
    end

    
    if( ~isfield( usr_par, 'use_mex') )
        usr_par.use_mex = 'no';
    end
    
    
    if( strcmp( usr_par.type, 'structure' ) )
        if( isfield( usr_par, 'structure_inversion') )
            if( ~isfield( usr_par.structure_inversion, 'rho') )
                [~,~,nx,nz,~,~,~,model_type] = input_parameters();
                [~,rho] = define_material_parameters( nx, nz, model_type);                
                usr_par.structure_inversion.rho = rho;
            end
            
            if( ~isfield( usr_par.structure_inversion, 'v0') )
                usr_par.structure_inversion.v0 = 4000;
            end
        else
            [~,~,nx,nz,~,~,~,model_type] = input_parameters();
            [~,rho] = define_material_parameters( nx, nz, model_type);
            usr_par.structure_inversion.rho = rho;
            
            usr_par.structure_inversion.v0 = 4000;
        end
    end

    
    if( ~isfield( usr_par, 'measurement') )
        usr_par.measurement = 'waveform_difference';
    end

    
    if( ~isfield( usr_par, 'veldis') )
        usr_par.veldis = 'dis';
    end

    
    if( isfield(usr_par,'filter') )
        if( ~isfield( usr_par.filter, 'apply_filter') )
            usr_par.filter.apply_filter = 'no';
        end
    else
        usr_par.filter.apply_filter = 'no';
    end

    
    if( ~isfield( usr_par, 'network') )
        usr_par.network = load('../output/interferometry/array_1_ref.mat');
    end

    
    if( ~isfield( usr_par, 'data') )
        usr_par.data = load('../output/interferometry/data_1_ref_0.mat');
    end
    
    
    if( isfield(usr_par,'kernel') )
        if( ~isfield( usr_par.kernel, 'percentile') )
            usr_par.kernel.percentile = 0;
        end
        
        if( ~isfield( usr_par.kernel, 'imfilter') )
            usr_par.kernel.imfilter = fspecial('gaussian',[60 60], 30);
        end
    else
        usr_par.kernel.percentile = 0;                
        usr_par.kernel.imfilter = fspecial('gaussian',[60 60], 30);
    end
    
    
    if( isfield(usr_par,'regularization') )
        if( ~isfield( usr_par.regularization, 'alpha') )
            usr_par.regularization.alpha = 0.0;
        end
        
        if( ~isfield( usr_par.regularization, 'weighting') )
            usr_par.regularization.weighting = 1.0;
        end
    else
        usr_par.regularization.alpha = 0.0;
        usr_par.regularization.weighting = 1.0;
    end
    
    
    if( isfield(usr_par,'debug') )
        if( ~isfield( usr_par.debug, 'df') )
            usr_par.debug.switch = 'no';
        end
        
        if( ~isfield( usr_par.debug, 'df') )
            usr_par.debug.df = 0;
        end
    else
        usr_par.debug.switch = 'no';
        usr_par.debug.df = 0;
    end


end
