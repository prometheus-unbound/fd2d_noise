
%==========================================================================
% loading function for parfor-loops (only for 1 variable)
%
% [x] = parload( filename )
%
% input:
%--------
% filename: name for file
%
% output:
%--------
% x: variable contained in file
%
%==========================================================================


function [x] = parload(filename)


    foo = load(filename);
    whichVariables = fieldnames(foo);

    if (numel(whichVariables) == 1)
        x = foo.(whichVariables{1});
    else
        error(['Problem with ', filename, '!']);
    end
    

end
