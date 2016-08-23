
%==========================================================================
% check if all necessary fields are set
%
% [usr_par] = usr_par_init_default_parameters(usr_par)
%
%==========================================================================


function [usr_par] = usr_par_init_default_parameters(usr_par)


    if (~isfield(usr_par, 'type'))
        usr_par.type = 'source';
    end


    if (isfield(usr_par, 'measurement'))
        if (~isfield(usr_par.measurement, 'type'))
            usr_par.measurement.type = 'waveform_difference';
        end
    else
        usr_par.measurement.type = 'waveform_difference';
    end


    if (isfield(usr_par, 'measurement'))
        if (~isfield(usr_par.measurement, 'mode'))
            usr_par.measurement.mode = 'auto';
        end
    else
        usr_par.measurement.mode = 'auto';
    end

    
    if (~isfield(usr_par, 'verbose'))
        usr_par.verbose = true;
    end
    

    if (isfield(usr_par, 'smoothing'))
        if (~isfield(usr_par.smoothing, 'sigma'))
            usr_par.smoothing.sigma = [10 10];
        end
    else
        usr_par.smoothing.sigma = [10 10];
    end


    if (~isfield(usr_par, 'data_independent'))
        usr_par.data_independent = 'no';
    end


    if (~isfield(usr_par, 'network'))
        error('\nload array!\n')
    end


    if (~isfield(usr_par, 'data'))
        if (strcmp(usr_par.data_independent, 'no'))
            error('\nload data!\n')
        end
    end


end

