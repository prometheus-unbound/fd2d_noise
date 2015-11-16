
function [j, g, c] = get_obj_grad(m, usr_par)
    
    [j, g, c] = eval_objective_and_gradient(m, 'yo', usr_par);

end


