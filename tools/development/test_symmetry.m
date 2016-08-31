
clear all
clc

nx = 20;
nz = 20;

v = rand(nx,nz);
w = rand(nx,nz);

x = linspace( 1, 500, nx );
z = linspace( 1, 500, nz );
[X, Z] = meshgrid(x, z);


test1 = norm(v' * ( gaussblur2d(w,x,z,[200 200]) ));
test2 = norm(w' * ( gaussblur2d(v,x,z,[200 200]) ));


mine = fspecial('gaussian',[10, 10], 30);
test3 = norm(v' * ( imfilter(w,mine,'circular') ));
test4 = norm(w' * ( imfilter(v,mine,'circular') ));