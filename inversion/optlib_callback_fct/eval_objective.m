
function [j, c] = eval_objective(m, ModRandString, usr_par)
    
    [j, ~, c] = eval_objective_and_gradient(m, ModRandString, usr_par);

end


