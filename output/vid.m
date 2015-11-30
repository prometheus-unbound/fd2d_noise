

for n = 1:size(G_2,3)
    
    mesh(G_2(:,:,end-n+1))
    drawnow
    pause(0.01)
    
end