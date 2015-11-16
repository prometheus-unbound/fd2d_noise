
function [Lx,Lz,nx,nz,dt,nt,v,rec_x,rec_z] = california()


%% select year 2005 for USArray and only westcoast
tstart = datetime(2005,1,1);
tend = datetime(2005,12,31);

lonlim = [-130 -110];
latlim = [30 45];

import_usarray
indices = find( ~isnat(starttime) & ~isnat(endtime) & starttime<tstart & endtime>tend & usalat>latlim(1) & usalat<latlim(2) & usalon>lonlim(1) & usalon<lonlim(2));
usalat_sel = usalat(indices);
usalon_sel = usalon(indices);


% read shapefile with land areas
land = shaperead('landareas', 'UseGeoCoords', true, 'BoundingBox', [lonlim', latlim']);

% map land onto regular grid
density = 100;
[field, R] = vec2mtx([land.Lat], [land.Lon], density, 'filled');

% map stations on regular grid
field_r = imbedm(usalat_sel, usalon_sel, 10, field, R);

% define X and Y grids
[X,Z] = meshgrid( linspace(min(land.Lon),max(land.Lon),size(field,2)) , linspace(min(land.Lat),max(land.Lat),size(field,1)) );

% remove great salt lake
[salt_lake_x, salt_lake_y] = find( X<-110 & X>-115 & Z>40 & Z<42 );
field(salt_lake_x,salt_lake_y) = 0.0;


%% select part of map
% Shapiro setup
% x = [-125 -113.5];
% z = [32 43];

x = [-126 -112.5];
z = [28 44];


% approximate factor to convert deg to m
x_deg2m = 40000/360 * 10^3 * sin( (z(1)+z(2))/2/180*pi );
y_deg2m = 40000/360 * 10^3;

% select specified area
[~, area_x] = find( X(1,:)>=x(1) & X(1,:)<=x(2) );
[area_z, ~] = find( Z(:,1)>=z(1) & Z(:,1)<=z(2) );

area_x = area_x(1,1:2:end);
area_z = area_z(1:1:end,1);



% convert deg to m
X_area = X(area_z, area_x) * x_deg2m;
Z_area = Z(area_z, area_x) * y_deg2m;


% set (0,0) to southwestern point
X_area = X_area - min(min(X_area));
Z_area = Z_area - min(min(Z_area));
area = field(area_z, area_x)';
area_r = field_r(area_z, area_x)';

% get x- and z-coordinates of receivers
[i] = find( area_r' == 10 );
rec_x = X_area(i);
rec_z = Z_area(i);


%% plot selected area
figure
pcolor(X_area,Z_area,area_r')
shading interp
axis image
view([0 90])

v = (area<=1)*3000 + (area==2)*3500; 


%% output variables
Lx = max(max(X_area));
Lz = max(max(Z_area));

nx = size(X_area,2);
nz = size(Z_area,1);

dt = 0.25 * min(Lx/(nx-1),min(Lz/(nz-1))) / max(max( v ));

npl = 15; 
fprintf( 'dx_max = %f\n', min(min(v)) / ( (npl-1) * 0.18 ) );
fprintf( 'freq_max = %f\n', min(min(v)) / ( (npl-1) * max( Lx/(nx-1), Lz/(nz-1) ) ) );
fprintf( 'lamda_min = %f\n\n', (npl-1) * max(Lx/(nx-1),Lz/(nz-1)) );

index = 1;
for i = 1:(size(rec_x,1)-1)
    for j = (i+1):size(rec_x,1)
        dist(index,1) = sqrt( ( rec_x(i) - rec_x(j) ).^2 + ( rec_z(i) - rec_z(j) ).^2 );
        index = index + 1;
    end
end

nt = max(dist) / min(min(v)) / dt + 40/dt;


end

