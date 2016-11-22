function [ j, g ] = regularization( j, g, m, usr_par )

    nx = usr_par.config.nx;
    nz = usr_par.config.nz;
    

    % misfit with regularization
    j = j + usr_par.regularization.alpha/2 * ...
        sum( repmat( usr_par.regularization.weighting, ( numel(m)-nx*nz ) / numel(usr_par.regularization.weighting), 1 ) .* ...
        ( m( 1 : end - nx*nz, 1 ) - usr_par.m0( 1 : end - nx*nz, 1 ) ).^2 );

    j = j + usr_par.regularization.beta/2 * ...
        sum( usr_par.regularization.weighting .* ...
        ( m( end - nx*nz + 1 : end, 1 ) - usr_par.m0( end - nx*nz + 1  : end, 1 ) ).^2 );
    
    
    % gradient with regularization
    if( ~ isscalar(g) )
        g( 1 : end - nx*nz, 1 ) = g( 1 : end - nx*nz, 1 ) + ...
            usr_par.regularization.alpha * repmat( usr_par.regularization.weighting, ( numel(m)-nx*nz ) / numel(usr_par.regularization.weighting), 1 ) .* ...
            ( m( 1 : end - nx*nz, 1 ) - usr_par.m0( 1 : end - nx*nz, 1 ) );
        
        g( end - nx*nz + 1 : end, 1 ) = g( end - nx*nz + 1 : end, 1 ) + ...
            usr_par.regularization.beta * usr_par.regularization.weighting .* ...
            ( m( end - nx*nz + 1 : end, 1 ) - usr_par.m0( end - nx*nz + 1 : end, 1 ) );
    end
    

end

