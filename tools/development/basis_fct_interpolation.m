
clear all
% close all
clc


n_basis_fct = 15;


%% calculate smooth spectrum
n_noise_sources = 1;
f_peak = [1/11, 1/7];
bandwidth = [0.035, 0.025];

strength = [0.5, 1];
[f_sample, n_sample] = input_interferometry();
spectrum = zeros(n_sample,n_noise_sources);
for ns = 1:n_noise_sources
    spectrum(:,ns) = strength(ns) * exp( -(abs(f_sample)-f_peak(ns)).^2 / bandwidth(ns)^2 );
end


%% compute basis function representation
[int_limits] = integration_limits(n_sample,n_basis_fct);
freq_basic = zeros(n_basis_fct,1);
spectrum_basic = zeros(n_basis_fct,1);
spectrum_box = zeros(n_sample,1);

for ib = 1:n_basis_fct
    
    indices = int_limits(ib,1) : int_limits(ib,2);
    ni = length(indices);
    
    freq_basic(ib,1) = f_sample( indices(1) + round( ni/2 ) - 1 );
    
    for k = indices
        for ns = 1:n_noise_sources
            spectrum_basic(ib,1) = spectrum_basic(ib,1) + spectrum(k,ns) / ni;
        end
    end
    
    for k = indices
        spectrum_box(k,1) = spectrum_box(k,1) + spectrum_basic(ib,1);
    end
    
end


%% reconstruct smooth version from basis function version
spectrum_reconst = zeros(n_sample,1);

for k = 1:n_sample
    
    ib = find( k >= int_limits(:,1) & k <= int_limits(:,2) );    
    indices = int_limits(ib,1) : int_limits(ib,2);
    ni = length(indices);
    
    current = spectrum_basic( ib );
    
    if( ib == 1 )
        below = 0;
        above = spectrum_basic( 2 );
    elseif( ib == n_basis_fct )
        below = spectrum_basic( n_basis_fct - 1 );
        above = 0;
    else
        below = spectrum_basic( ib - 1 );
        above = spectrum_basic( ib + 1 ); 
    end
    
    i = find( indices == k );
    if( i < (ni+1)/2 )        
        spectrum_reconst(k,1) = current + ((ni+1)/2-i) * ( below - current ) / ni;
    elseif( i == (ni+1)/2 )
        spectrum_reconst(k,1) = current;        
    elseif( i > (ni+1)/2 )
        spectrum_reconst(k,1) = current + (i-(ni+1)/2) * ( above - current ) / ni;
    end
    
end



matrix = zeros( n_sample, n_basis_fct );
for k = 1:n_sample
    
    ib = find( k >= int_limits(:,1) & k <= int_limits(:,2) );    
    indices = int_limits(ib,1) : int_limits(ib,2);
    ni = length(indices);
    
    current = ib;
    
    if( ib == 1 )
        below = [];
        above = 2;
    elseif( ib == n_basis_fct )
        below = ib-1;
        above = [];
    else
        below = ib - 1;
        above = ib + 1; 
    end
    
    i = find( indices == k );
    if( i < (ni)/2 )        
        
        matrix( k, current ) = 1 - ((ni)/2-i)/ni;
        if( ~isempty(below) )
            matrix( k, below ) = ((ni)/2-i)/ni;
        end
        
    elseif( i == (ni)/2 )
        
        matrix( k, current ) = 1;
        
    elseif( i > (ni)/2 )
        
        matrix( k, current ) = 1 - (i-(ni)/2)/ni;
        if( ~isempty(above) )
            matrix( k, above ) = (i-(ni)/2)/ni;
        end
        
    end
    
end





%% plotting
figure(1)
clf
hold on
plot( f_sample, sum(spectrum,2), 'k' )
plot( f_sample, spectrum_box, 'bx' )
% plot( f_sample, spectrum_reconst, 'ro' )
plot( freq_basic, spectrum_basic, 'g+' )
plot( f_sample, matrix * spectrum_basic, 'ro')





