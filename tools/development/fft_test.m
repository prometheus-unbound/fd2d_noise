
clear all
clc


%% load array and data
load ../output/interferometry/array_16_ref.mat
data = load( '../output/interferometry/data_16_ref_91_2h2g_homog.mat' );
initial = load( '../output/interferometry/data_16_ref_91_2h_homog.mat' );


%% select one correlation
n_ref = size( array, 1 );
n_rec = n_ref-1;
t = data.t;

Nsamps = length(t);
dt = t(2)-t(1);
Fs = 1/dt;

i_ref = 3;
i_rec = 10;
i_data = (i_ref-1) * n_rec + i_rec;

src = ref_stat( i_ref, : );
rec = array( find( ~ismember(array, src, 'rows') ) ,:);


%% set up window
% left = t(1);
left = 0;
right = t(end);

% distance = sqrt( (src(1,1) - rec(i_rec,1)).^2 + (src(1,2) - rec(i_rec,2)).^2 );
% left = distance/4000.0 - 27.0;
% right = distance/4000.0 + 27.0;
% if( left < 0 )
%     index = find( t==0 );
%     left = t(index+1);
% end
% if( right > t(end) )
%     right = t(end);
% end

% win = ones( 1, length(t) );
win = get_window( t, left, right, 'cos_taper' );
% win = get_window( t, left, right, 'hann' );



%% compute adjoint source time function
y = data.c_data(i_data,:);

[~, adstf] = make_adjoint_sources( initial.c_data(i_data,:), y, 0*y, t, 'dis', 'waveform_difference', src, rec(i_rec,:), '1st' );
% [~, adstf] = make_adjoint_sources( initial.c_data(i_data,:), y, 0*y, t, 'dis', 'log_amplitude_ratio', src, rec(i_rec,:), '1st' );
adstf = fliplr( adstf );


% y = 0 * data.c_data(i_data,:);
% y(1,1:Nsamps-1) = diff( data.c_data(i_data,:) ) / dt;
y = win .* y;
adstf = win .* adstf;





% indices = find( t >= 160 & t <= 180 );
% y( indices ) = y( indices ) / 20;



% Fourier Transform
f = Fs*(0:Nsamps/2-1)/Nsamps;

data_fft = abs(fft(y));
data_fft = data_fft( 1:Nsamps/2 );

adstf_fft = abs(fft(adstf));
adstf_fft = adstf_fft( 1:Nsamps/2 );



%% plotten
figure(1)
clf
hold on
% plot(t, y, 'b')
plot(t, adstf, 'r')
% plot(t, win, 'k')
% plot(t, initial.c_data(i_data,:), 'c')
% plot(t, win.^2 .* ( initial.c_data(i_data,:) - y ) * dt, 'y')
xlabel('Time (s)')
ylabel('Amplitude')


figure(2)
clf
hold on
plot(f, data_fft / max(data_fft), 'b')
plot(f, adstf_fft / max(adstf_fft), 'r')
xlim([0 0.5])
xlabel('Frequency (Hz)')
ylabel('Amplitude')



