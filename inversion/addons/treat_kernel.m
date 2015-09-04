
function [K] = treat_kernel(K_raw,type,myfilter,percentile)
    

    [~,~,~,~,~,~,~,~,~,n_basis_fct] = input_parameters();
    
    
    %% clip kernel if desired
    if( percentile ~= 0 )
        
        if( n_basis_fct == 0 || strcmp(type,'structure') )
            n_basis_fct = 1;
        end
        
        K = K_raw;
        clip_value = zeros(1,n_basis_fct);
        for k = 1:n_basis_fct
        
            clip_value(k) = prctile( reshape( abs(K_raw(:,:,k)), [], 1), percentile );
            [i,j] = find( abs(K_raw(:,:,k)) >= clip_value(k) );
            
            K(i,j,k) = clip_value(k);
        
        end
        
    end
    
    
    %% apply filter according to myfilter 
    if( myfilter ~= 0 )
        
        K = imfilter( K_raw, myfilter, 'replicate' );
        
    end
    

end