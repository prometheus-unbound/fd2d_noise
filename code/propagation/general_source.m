
function [dist_general] = general_source( noise_spectrum, noise_source_distribution )
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % could be done in a nicer way
    %
    % keep number of degrees of freedom - mapping function from small setup
    % to map for each frequency
    % still thinking about it... ;)
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
            ni = length(indices);
            
            % combine to frequency bands - leads to the same source map for several frequencies
            % normalize the summation over frequencies by number of frequencies involved
            
            % % slow, but nicer version
            % for k = indices
            %     for ns = 1:n_noise_sources
            %         dist_general(:,:,indices) = dist_general(:,:,indices) + ...
            %             repmat( noise_spectrum(k,ns) * noise_source_distribution(:,:,ns) / ni, 1, 1, ni ) ;
            %     end
            % end
            
            % faster version
            dist_basis = zeros(nx,nz);
            for k = indices
                for ns = 1:n_noise_sources
                    dist_basis = dist_basis + noise_spectrum(k,ns) * noise_source_distribution(:,:,ns);
                end
            end
            
            for k = indices
                dist_general(:,:,k) = dist_basis/ni;
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