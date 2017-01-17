
clear all


n_eigenvectors = 100;


load('../output/interferometry/data_16_ref_0_gaussian_random_0.07_0.8e10_nosmooth.mat')
c_data_all = zeros( 240*n_eigenvectors, 2599 );
c_data_new = zeros( 240, 2599 );
weights = zeros(100,1);


for i = 1:n_eigenvectors
   
    c_data_all( 1+(i-1)*240 : i*240, : ) = c_data;
    
    if( mod(i,2) == 1 )
        weights(i,1) = -1;
    else
        weights(i,1) = 1;
    end
    
end


tic
for i = 1:n_eigenvectors
   
    c_data_new = c_data_new + weights(i,1) * c_data_all( 1+(i-1)*240 : i*240, : );
    
end
toc