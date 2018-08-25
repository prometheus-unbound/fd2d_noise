
clear all
% close all

% load ~/Desktop/data_1_ref_0_homogeneous.mat
% load ~/Desktop/data_1_ref_0_point.mat
% load ~/Desktop/data_1_ref_0_line3.mat
% load ~/Desktop/data_1_ref_0_line3_m.mat
load ~/Desktop/data_1_ref_0_line3_m3.mat
% load ~/Desktop/data_1_ref_0_line4.mat

correlation = filter_correlations(c_data,t,0.05,0.18,4);


fig1 = figure(1);
set(fig1, 'units', 'normalized', 'position', [0.1 0.6 0.5 0.13])

plot(t,  correlation/max(abs(correlation)), 'k', 'linewidth', 2)
hold on

font_size = 18;
title_size = 24;
maker_size = 10;
set(gca, 'FontSize', font_size);
     
xlabels = [-150 -75 0 75 150];
ylabels = [-1 0 1];
set(gca, 'XTick', xlabels);
set(gca, 'YTick', ylabels);
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);

xlabel('t [s]')
% ylabel('normalized')
% title('homogeneous distribution', 'FontSize', title_size)
title('noise anomalies 1-7', 'FontSize', title_size)
% title('noise anomaly 1', 'FontSize', title_size)


xlim([-160 160])
ylim([-1 1])

ax = gca;
ax.LineWidth = 2;

orient landscape