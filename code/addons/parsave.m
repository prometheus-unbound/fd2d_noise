
function parsave(fname, x, y)
    
    if( isempty(y) )
        save(fname, 'x')
    else
        save(fname, 'x', 'y')
    end
    
end
