
load ../output/interferometry/data_16_ref_0_2h2g_nover.mat

Nsamps = length(t);

y = c_data(10,:);
y = filter_correlations( y, t, 1/7-0.01, 1/7+0.01 );
y = filter_correlations( y, t, 1/7-0.01, 1/7+0.01 );
y = filter_correlations( y, t, 1/7-0.01, 1/7+0.01 );
y = filter_correlations( y, t, 1/7-0.01, 1/7+0.01 );

dt = t(2)-t(1);
Fs = 1/dt;

%Do Fourier Transform
y_fft = abs(fft(y));            %Retain Magnitude
y_fft = y_fft(1:Nsamps/2);      %Discard Half of Points
f = Fs*(0:Nsamps/2-1)/Nsamps;   %Prepare freq data for plot

%Plot Sound File in Time Domain
figure
plot(t, y)
xlabel('Time (s)')
ylabel('Amplitude')
title('Tuning Fork A4 in Time Domain')

%Plot Sound File in Frequency Domain
figure
plot(f, y_fft)
% xlim([0 1000])
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Frequency Response of Tuning Fork A4')