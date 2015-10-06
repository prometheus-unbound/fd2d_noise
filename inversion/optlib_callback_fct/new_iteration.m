function [usr_par] = new_iteration( it, model, usr_par) 
% NEW_ITERATION This auxiliary function is called before a new iteration is started.
% You might use it for writing output files, plotting, ...
%
%
% Input: 
% m : current model
% it : current iteration number
% usr_par : auxiliary user defined parameters (optional)
%
% Output:
% usr_par : auxiliary user defined parameters (optional)


if(it>0)
    model.normg_frac = norm(model.gradient) / usr_par.model(1).normg;
else
    model.normg_frac = 1.0;
end

fprintf('\n%i\t misfit = %e\t norm(g)/norm(g0)=%e\n', it, model.objective, model.normg_frac)

usr_par.model(it+1) = orderfields( model );

end
