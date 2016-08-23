
function [g] = eval_grad_objective(m, ModRandString, usr_par)
    
    [~, g] = eval_objective_and_gradient(m, ModRandString, usr_par);

end


