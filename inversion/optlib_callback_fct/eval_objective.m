
function [j] = eval_objective(m, ModRandString, usr_par)
    
    [j, ~] = eval_objective_and_gradient(m, ModRandString, usr_par);

end


