
clear all
% close all
clc


green = load('~/Desktop/green_1_ref_1.mat');
% green = load('../output/interferometry/data_1_ref_0.mat');
% corr = load('~/Desktop/correlation_1_ref_2_heto2.mat');
load ~/Desktop/array_1_ref.mat
% load ../output/interferometry/array_1_ref.mat


dt = green.t(2)-green.t(1);
green.t = 0:dt:(size(green.c_data,2)-1)*dt;
% plot_recordings(green.c_data,green.t,'dis','b',false)
% plot_recordings(corr.c_data,corr.t,'vel','r',true)

     
for k = 1:size(green.c_data,1)
    green.c_data(k,1:end-1) = diff( green.c_data(k,:) ) / dt;
    green.c_data(k,end) = 0.0;
end

green.c_data = filter_correlations( green.c_data, green.t, 0.01, 0.1, 4 );  
% plot_recordings(green.c_data,green.t,'dis','b',true)

% for k = 1:size(corr.c_data,1)
%     corr.c_data(k,1:end-1) = diff( corr.c_data(k,:) ) / dt;
%     corr.c_data(k,end) = 0.0;
% end


% i_zero = find( corr.t==0 );
% amp_corr = zeros( size(corr.c_data,1)-1, 1 );
amp_green = zeros( size(green.c_data,1)-1, 1 );
distance = zeros( size(green.c_data,1)-1, 1 );

% 1
src_x = 500834.7245409015;
rec_x = [601001.6694490818,   701168.6143572620,   801335.5592654423,   901502.5041736227,   998330.5509181969,  1098497.4958263773,  1198664.4407345578,  1298831.3856427381,  1398998.3305509184,  1499165.2754590986];

% 2
% src_x = 500625.7822277847;
% rec_x = [600750.9386733416,   700876.0951188987,   801001.2515644556,   901126.4080100126,   998748.4355444306,  1098873.5919899875,  1198998.7484355443,  1299123.9048811013,  1399249.0613266584,  1499374.2177722151];

% 3
% src_x = 500500.5005005005;
% rec_x = [600600.6006006006,   700700.7007007007,   800800.8008008008,   900900.9009009008,   998998.9989989990,  1099099.0990990992,  1199199.1991991992,  1299299.2992992993,  1399399.3993993993,  1499499.4994994996];

% 4 5
% src_x = 501253.1328320802;
% rec_x = [601503.7593984962,   701754.3859649124,   802005.0125313284,   902255.6390977445,   997493.7343358396,  1097744.3609022554,  1197994.9874686715,  1298245.6140350876,  1398496.2406015038,  1498746.8671679199];

for i = 1:length(rec_x)
    
    amp_green( i ) = max( abs(green.c_data( i, : )) );
%     amp_corr( i ) = max( abs(corr.c_data( i+1, i_zero:end )) );
    
    distance( i ) = rec_x( i ) - src_x( 1 );
    
end


figure(1)
% semilogy( distance, amp_corr/max(amp_corr), 'rx--' )

semilogy( distance, amp_green / max(amp_green), 'r+--' )
% semilogy( distance, amp_green / mean(amp_green), 'r+--' )
hold on
test = 1 ./ (distance).^(1/2);
semilogy( distance, test / max( test ), 'ko--' )
% semilogy( distance, test/mean(test), 'ko--' )

% legend('correlation', 'green', 'analytical')
legend('green', 'analytical')


% figure
% plot( (distance).^(1/2) .* amp_green / max( (distance).^(1/2) .* amp_green )  )


