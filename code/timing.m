
for i = 1:n_ref
        
    src = ref_stat(i,:);
    rec = array( find(~ismember(array,src,'rows') ) , :);
    
    fprintf( 'ref %i: calculate Green function\n', i )    
    G_fft = run_forward1_green( structure, src, 0 );
    
end