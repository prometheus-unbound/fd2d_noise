
function [dist_general] = general_source( noise_spectrum, noise_source_distribution )
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % could be done in a nicer way
    %
    % keep number of degrees of freedom - mapping function from small setup
    % to map for each frequency
    % don't want to keep noise spectrum separate, thinking about it... ;)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    [~,~,nx,nz,~,~,~,~,~,n_basis_fct] = input_parameters();
    [~,n_sample] = input_interferometry();
    n_noise_sources = size( noise_spectrum, 2 );
   
    
    % source map for each frequency - done as in run_forward
    dist_general = zeros(nx,nz,n_sample);
    
    if( n_basis_fct ~= 0 )
        
        int_limits = integration_limits(n_sample,n_basis_fct);
        
        for ib = 1:n_basis_fct
            
            indices = int_limits(ib,1) : int_limits(ib,2);
            
            % combine to frequency bands - leads to the same source map for several frequencies
            % normalize the summation over frequencies by number of frequencies involved
            for k = indices
                for ns = 1:n_noise_sources
                    dist_general(:,:,indices) = dist_general(:,:,indices) + ...
                        repmat( noise_spectrum(k,ns) * noise_source_distribution(:,:,ns), 1, 1, length(indices) ) / length(indices);
                end
            end
            
        end
        
    else
        
        % combine noise spectrum and noise source distribution
        for k = 1:n_sample
            for ns = 1:n_noise_sources
            
                dist_general(:,:,k) = dist_general(:,:,k) + noise_spectrum(k,ns) * noise_source_distribution(:,:,ns);
                
            end
        end
        
    end
    
    
end