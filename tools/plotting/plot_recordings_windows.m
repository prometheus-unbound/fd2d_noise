
function plot_recordings_windows(u,t,veldis,color,normalize,left,right)

%- initialisations and parameters -----------------------------------------

spacing = 2;
a = 0;


%- convert to velocity if wanted ------------------------------------------

if strcmp(veldis,'vel')

    nt = length(t);
    v = zeros(size(u,1),nt-1);
    
    for k=1:size(u,1)
        v(k,:) = diff(u(k,:)) / (t(2)-t(1));
    end
    
    t = t(1:nt-1);
    u = v;

end

%- plot recordings with ascending distance from the first source ----------

% figure
set(gca,'FontSize',20)
hold on

for k=1:size(u,1)
    
    if(normalize == true)
        m = max(abs(u(k,:)));
    else
        m = 1;
    end
    
    plot(t,spacing*(k-1+a)+u(k,:)/m,color,'LineWidth',1);
    
    plot([left(k,1) left(k,1)], [spacing*(k-2+a)+0.5 spacing*(k+a)-0.5],'b--')
    plot([right(k,1) right(k,1)], [spacing*(k-2+a)+0.5 spacing*(k+a)-0.5],'b--')
    
    tmp = left;
    left = -right;
    right = -tmp;
    clear tmp;
    
    plot([left(k,1) left(k,1)], [spacing*(k-2+a)+0.5 spacing*(k+a)-0.5],'b--')
    plot([right(k,1) right(k,1)], [spacing*(k-2+a)+0.5 spacing*(k+a)-0.5],'b--')
    
end

xlabel('time [s]','FontSize',20);
ylabel('normalised traces','FontSize',20);
% set(gca, 'YTick', []);
