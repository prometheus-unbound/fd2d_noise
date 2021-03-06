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


%- compute ratio of gradients ---------------------------------------------
if(it>0)
    model.normg_frac = norm(model.gradient) / usr_par.model(1).normg;
else
    model.normg_frac = 1.0;
end


%- some output - a bit puristic at the moment -----------------------------
fprintf('\n%i\t misfit = %e\t norm(g)/norm(g0)=%e\n', ...
    it, model.objective, model.normg_frac)


%- save model structure in usr_par structure ------------------------------
usr_par.model(it+1) = orderfields( model );


%- save model -------------------------------------------------------------
% model.imfilter.source = usr_par.kernel.imfilter.source;
% model.imfilter.structure = usr_par.kernel.imfilter.structure;

model.sigma.source = usr_par.kernel.sigma.source;
model.sigma.structure = usr_par.kernel.sigma.structure;
model.config.n_basis_fct = usr_par.config.n_basis_fct;
model.config.nx = usr_par.config.nx;
model.config.nz = usr_par.config.nz;

if( strcmp( usr_par.ring.switch, 'yes' ) )
   model.ring = usr_par.ring; 
end

save( sprintf('../inversion/models/model_%i.mat',it), 'model' )


end
