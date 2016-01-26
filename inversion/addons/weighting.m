function [ w ] = weighting( nx, nz, n_basis_fct )

    if( n_basis_fct == 0 )
        n_basis_fct = 1;
    end

    % set weighting within domain to 4
    w = 4 * ones(nx, nz, n_basis_fct + 1);
    
    % set weighting at boundaries to 2
    w (1,:,:) = 2;
    w (end,:,:) = 2;
    w (:,1,:) = 2;
    w (:,end,:) = 2;
    
    % set weighting in corners to 1
    w(1,1,:) = 1;
    w(1,end,:) = 1;
    w(end,1,:) = 1;
    w(end,end,:) = 1;
    
    % reshape to vector
    w = reshape( w, [], 1 );

end

