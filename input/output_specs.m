 

%==========================================================================
% print output concerning specification and timing
%==========================================================================

verbose = 'no';


%==========================================================================
% make plots every nth timestep
%==========================================================================

make_plots = 'no';
plot_nt = 100;


%==========================================================================
% summarize configuration
%==========================================================================

if( strcmp(verbose,'yes') )
    
    v_clf = 4000;
    fprintf('\ndx = %f\n',Lx/(nx-1));
    fprintf('dz = %f\n',Lz/(nz-1));
    fprintf('dt = c * min(dx,dz)/%5.1f = c * %f\n',v_clf,min(Lx/(nx-1),Lz/(nz-1))/v_clf);
    fprintf('chosen c = %f\n\n',dt/(min(Lx/(nx-1),Lz/(nz-1))/v_clf));
    fprintf('freq_max = %f\n',v_clf/(15*max(Lx/(nx-1),Lz/(nz-1))));
    fprintf('lamda_min = %f\n\n',15*max(Lx/(nx-1),Lz/(nz-1)));
    
    fprintf('freq_min = %f\n', v_clf/(sqrt((rec_x-src_x).^2 + (rec_z-src_z).^2)/3) );
    fprintf('lamda_max = %f\n', sqrt((rec_x-src_x).^2 + (rec_z-src_z).^2)/3 );
    
end


