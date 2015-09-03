
function [c_data] = filter_correlations(c_data,t,f_min,f_max)
    
    for i = 1:size(c_data,1)
        
        c_data(i,:)  = fliplr( butterworth_lp( fliplr( c_data(i,:) ), t, 5, f_max, 'silent') );
        c_data(i,:)  = butterworth_lp( c_data(i,:), t, 5, f_max, 'silent');

        c_data(i,:)  = fliplr( butterworth_hp( fliplr( c_data(i,:) ), t, 3, f_min, 'silent') );
        c_data(i,:)  = butterworth_hp( c_data(i,:), t, 3, f_min, 'silent');
        
    end

end