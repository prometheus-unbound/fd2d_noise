
%==========================================================================
% saving function for parfor-loops (only for one variable)
%
% parsave( filename, x )
%
% input:
%--------
% filename: name for file
% x: variable
%
%==========================================================================


function parsave(filename, x)

    save(filename, 'x')

end

