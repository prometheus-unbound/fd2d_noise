function [ w ] = weighting( nx, nz )

    % set weighting within domain to 4
    w = 4 * ones(nx, nz);
    
    % set weighting at boundaries to 2
    w (1,:) = 2;
    w (end,:) = 2;
    w (:,1) = 2;
    w (:,end) = 2;
    
    % set weighting in corners to 1
    w(1,1) = 1;
    w(1,end) = 1;
    w(end,1) = 1;
    w(end,end) = 1;
    
    % reshape to vector
    w = reshape( w, [], 1 );

end

