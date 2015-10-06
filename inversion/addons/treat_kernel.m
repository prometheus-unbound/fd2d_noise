
function [K] = treat_kernel( K_raw, usr_par )
    

    [~,~,~,~,~,~,~,~,~,n_basis_fct] = input_parameters();
    
    
    % clip kernel if desired
    if( usr_par.kernel.percentile ~= 0 )
        
        if( n_basis_fct == 0 || strcmp( usr_par.type, 'structure') )
            n_basis_fct = 1;
        end
        
        clip_value = zeros(1,n_basis_fct);
        for k = 1:n_basis_fct
        
            clip_value(k) = prctile( reshape( abs(K_raw(:,:,k)), [], 1), percentile );
            [i,j] = find( abs(K_raw(:,:,k)) >= clip_value(k) );
            
            K(i,j,k) = clip_value(k);
        
        end
            
    
    % apply filter according to myfilter 
    elseif( usr_par.kernel.smoothing ~= 0 )
        
        K = imfilter( K_raw, usr_par.kernel.smoothing, 'replicate' );
        
        
    else
        
        K = K_raw;
        
    end
    

end