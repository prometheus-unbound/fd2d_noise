% plot_recordings(u,t,veldis,color,normalize)
%
% u: displacement recordings
% t: time axis
% veldis: 'dis' for displacement, 'vel' for velocity
% color
% normalize, true or false

function h = plot_recordings(u,t,veldis,color,normalize)

spacing = 2;
a = 0;

%- convert to velocity if wanted ------------------------------------------
nt = length(t);

if strcmp(veldis,'vel')

    v = zeros(size(u,1),nt-1);
    
    for k = 1:size(u,1)
        v(k,:) = diff(u(k,:)) / (t(2)-t(1));
    end
    
    t = t(1:nt-1);
    u = v;

end


%- filter displacement recordings -----------------------------------------
% for i_rec = 1:size(u, 1)
%     u(i_rec, :) = filter_correlations( u(i_rec, :), t, 0.02, 0.2, 1 );
% end
    

%- plot recordings with ---------------------------------------------------
set(gca,'FontSize',20)
hold on

for k=1:size(u,1)
    
    % normalization
    if( normalize == true )
        m = max(abs(u(k,:)));
    else
        m=1;
    end
    
    % plot recordings
    h = plot(t,spacing*(k+a)+u(k,:)/m,'color',color,'LineWidth',1);
       
end

xlabel('time [s]','FontSize',20);
ylabel('normalised traces','FontSize',20);
set(gca, 'YTick', []);
