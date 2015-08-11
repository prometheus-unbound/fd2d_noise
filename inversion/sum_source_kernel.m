
function [K_s] = sum_source_kernel(K_all)
    
    [~,~,nx,nz,~,~,~,~,~,n_basis_fct] = input_parameters();
    [~,n_sample] = input_interferometry();
    
    
    % same update for each frequency
    if( n_basis_fct == 0 )
        
        K_s = repmat( sum(K_all,3), 1, 1, n_sample );
        
        
    % updates in different frequency bands
    else
        
        K_s = zeros(nx,nz,n_sample);
        [int_limits] = integration_limits(n_sample,n_basis_fct);
        
        for ib = 1:n_basis_fct
            
            indices = int_limits(ib,1) : int_limits(ib,2);
            
            % possibilities in changing the way how K_basis is built
            % try "overlapping" frequency bands, i.e. K_basis includes also neighbours
            K_basis = sum( K_all(:,:,indices), 3 ) / length(indices);
            
            K_s(:,:,indices) = repmat( K_basis, 1 , 1, length(indices) );
            
        end
        
    end
    

end