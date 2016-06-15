
clear all

checker_size = 50;
checkerboard  = zeros(600,600);

myfilter = fspecial('gaussian',[70 70], 50);
mesh(myfilter)

prev_i = 0;
prev_j = 0;

k = 1;
l = 1;

for i = 1:size(checkerboard,1)
    for j = 1:size(checkerboard,2)
        
        if( i < 2*checker_size || j < 2*checker_size || i > 600-2*checker_size || j > 600-2*checker_size)
            continue
        end
        
        if( floor((i-1)/checker_size) > prev_i )
            k = -k;
        end
        
        if( floor((j-1)/checker_size) > prev_j )
            l = -l;
        end
        
        prev_i = floor((i-1)/checker_size);
        prev_j = floor((j-1)/checker_size);
        
        checkerboard(i,j) = k*l;
    end
    
    prev_j = 0;
    l = 1;
end

checkerboard_filtered = imfilter( checkerboard, myfilter, 'symmetric' );


% figure(1)
% clf
% pcolor(checkerboard)
% shading interp
% axis square

figure(2)
clf
mesh(checkerboard_filtered)
shading interp
axis square
% view([0 90])
