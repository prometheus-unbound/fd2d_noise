
function [int_limits] = integration_limits(n_sample,n_basis_fct)

    step_length = n_sample / n_basis_fct;
    
    int_limits = zeros(n_basis_fct,2);
    for i = 1:n_basis_fct
        
        int_limits(i,1) = floor( (i-1) * step_length + 1 );
        int_limits(i,2) = floor( i * step_length );
        
    end

end

