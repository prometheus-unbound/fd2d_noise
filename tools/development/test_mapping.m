
clear all
clc


nt = 40;
dt = 0.1;

t = 0:dt:(nt-1)*dt;
m = sinc( 2*t )';


figure(1)
hold on
plot( t, m, 'k' )



% set up mapping matrix
map = zeros( 2*nt-1, nt );

j = 0;
for i = 1:size( map, 1 )
    
    
    if( i < nt )
        map( i, end - j ) = 1;
        j = j + 1;
    elseif( i==nt )
        map( i, 1 ) = 1;
        j = 2;
    else
        map( i, j ) = 1;
        j = j + 1;
    end
    
end

m_parameter = map * m;
t2 = -(nt-1)*dt:dt:(nt-1)*dt;
plot( t2, m_parameter, 'r' )


m2 = map' * m_parameter;

% wie ist das? im grunde ist man ja an der gesamten sensitvität
% interessiert und nicht an der durchschnittlichen
% m2(2:end) = m2(2:end) / 2;

plot( t, m2, 'b' )
