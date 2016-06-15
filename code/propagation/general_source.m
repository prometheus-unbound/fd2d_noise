
function [dist_basis_fct] = general_source( noise_spectrum, noise_source_distribution, n_basis_fct )

    
    [~,~,nx,nz] = input_parameters();
    [~,n_sample] = input_interferometry();
    n_noise_sources = size( noise_spectrum, 2 );
   
        
    dist_basis_fct = zeros(nx,nz,n_basis_fct);
    int_limits = integration_limits(n_sample,n_basis_fct);
    
    for ib = 1:n_basis_fct
        
        indices = int_limits(ib,1) : int_limits(ib,2);
        ni = length(indices);
        
        for k = indices
            for ns = 1:n_noise_sources
                dist_basis_fct(:,:,ib) = dist_basis_fct(:,:,ib) + noise_spectrum(k,ns) * noise_source_distribution(:,:,ns) / ni;
            end
        end
        
    end
    

end