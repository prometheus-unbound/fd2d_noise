
function [K] = treat_kernel( K_raw, usr_par )
       
    
    % clip kernel if desired
    if( usr_par.kernel.percentile ~= 0 )
        
        if( usr_par.config.n_basis_fct == 0 || strcmp( usr_par.type, 'structure') )
            n_basis_fct = 1;
        end
        
        clip_value = zeros(1,n_basis_fct);
        for k = 1:n_basis_fct
        
            clip_value(k) = prctile( reshape( abs(K_raw(:,:,k)), [], 1), percentile );
            [i,j] = find( abs(K_raw(:,:,k)) >= clip_value(k) );
            
            K(i,j,k) = clip_value(k);
        
        end
            
        
    else
        
        K = K_raw;
        
    end
    

end