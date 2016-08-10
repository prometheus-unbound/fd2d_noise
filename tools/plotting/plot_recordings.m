% plot_recordings(u, t, veldis, color, normalize)
%
% u: displacement recordings
% t: time axis
% color
% normalize, true or false

function h = plot_recordings(u,t,color,normalize)

spacing = 2;
offset = 0;

%- convert to velocity ----------------------------------------------------
nt = length(t);
v = zeros( size(u,1), size(u,2), nt-1);

for i_ref = 1:size(u,1)
    for i_rec = 1:size(u,2)
        v(i_ref,i_rec,:) = diff(u(i_ref,i_rec,:)) / (t(2)-t(1));
    end
end

t = t(1:nt-1);
u = v;


%- plot recordings --------------------------------------------------------
set(gca,'FontSize',18)
hold on
k = 1;
for i_ref=1:size(u,1)
    
    for i_rec=1:size(u,2)
        
        % normalization
        if( normalize == true )
            m = max(max(abs(u(i_ref,:,:))));
        else
            m = 1;
        end
        
        % plot
        h = plot(t, spacing * k + squeeze( u(i_ref,i_rec,:) ) / m + offset, 'color', color, 'LineWidth', 1);
        k = k+1;
    end
       
end

xlabel('time [s]','FontSize',18);
ylabel('normalised traces','FontSize',18);
set(gca, 'YTick', []);