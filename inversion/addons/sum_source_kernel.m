
function [K_s] = sum_source_kernel(K_all)
    
    [~,~,nx,nz,~,~,~,~,~,n_basis_fct] = input_parameters();
    [~,n_sample] = input_interferometry();
    
    
    % same update for all frequencies
    if( n_basis_fct == 0 )
        
        K_s = sum(K_all,3) / n_sample;
        
        
    % updates in different frequency bands
    else
        
        K_s = zeros(nx,nz,n_basis_fct);
        [int_limits] = integration_limits(n_sample,n_basis_fct);
        
        for ib = 1:n_basis_fct
            
            indices = int_limits(ib,1) : int_limits(ib,2);
            K_s(:,:,ib) = sum( K_all(:,:,indices), 3 ) / length(indices);
            
        end
        
    end
    

end