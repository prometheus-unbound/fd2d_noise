
function [m_parameters] = map_m_to_parameters( m, usr_par )
% MAP_M_TO_PARAMETERS function to map the model variable m to physical
% parameters.
%
% Input: 
% m
% usr_par : auxiliary user defined parameters (optional)
%
% Output:
% m_parameters
%
% See also MAP_PARAMETERS_TO_M and MAP_GRADPARAMETERS_TO_GRADM.


m_parameters = reshape( m, usr_par.config.nx, usr_par.config.nz, [] );

[absbound] = init_absbound();
[ix,iz] = find( absbound == 1, 1, 'first' );


if( strcmp( usr_par.ring.switch, 'yes' ) )
    
    [Lx, Lz] = input_parameters();
    [X, Z] = define_computational_domain( Lx, Lz, usr_par.config.nx, usr_par.config.nz );
    
    R = ( ( X - usr_par.ring.x_center_ring ).^2 + ( Z - usr_par.ring.z_center_ring ).^2 ).^(1/2);
    
    m_parameters(:,:,1) = m_parameters(:,:,1) .*  exp( -abs( R - usr_par.ring.radius ).^2 / usr_par.ring.taper_strength ) .* double(R > (usr_par.ring.radius-usr_par.ring.thickness/2) & R < (usr_par.ring.radius+usr_par.ring.thickness/2) );
    m_parameters(ix:end-ix+1, iz:end-iz+1, end) = imfilter( m_parameters(ix:end-ix+1, iz:end-iz+1, end) , usr_par.kernel.imfilter.structure, 'circular' );
    
else
    
    % m_parameters(ix:end-ix+1, iz:end-iz+1, 1:end-1) = imfilter( m_parameters(ix:end-ix+1, iz:end-iz+1, 1:end-1) , usr_par.kernel.imfilter.source, 'circular' );
    % m_parameters(ix:end-ix+1, iz:end-iz+1, end) = imfilter( m_parameters(ix:end-ix+1, iz:end-iz+1, end) , usr_par.kernel.imfilter.structure, 'circular' );
    
    [Lx, Lz] = input_parameters();
    [~, ~, x, z] = define_computational_domain(Lx, Lz, usr_par.config.nx, usr_par.config.nz);
    
    for i = 1:size(m_parameters,3)-1
        m_parameters(ix:end-ix+1, iz:end-iz+1, i) = gaussblur2d( m_parameters(ix:end-ix+1, iz:end-iz+1, i), x(ix:end-ix+1), z(iz:end-iz+1), usr_par.kernel.sigma.source );
    end
    m_parameters(ix:end-ix+1, iz:end-iz+1, end) = gaussblur2d( m_parameters(ix:end-ix+1, iz:end-iz+1, end), x(ix:end-ix+1), z(iz:end-iz+1), usr_par.kernel.sigma.structure );
    
end


m_parameters(:,:,1:end-1) = usr_par.initial.ref_source + m_parameters(:,:,1:end-1);
m_parameters(:,:,end) = usr_par.initial.mu_0 * ( usr_par.initial.ref_structure + m_parameters(:,:,end) );


end