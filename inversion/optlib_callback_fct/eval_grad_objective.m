
function [g, c] = eval_grad_objective(m, ModRandString, usr_par)
    
    [~, g, c] = eval_objective_and_gradient(m, ModRandString, usr_par);

end


