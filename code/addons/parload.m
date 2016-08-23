
%==========================================================================
% loading function for parfor-loops (only for 2 variables)
%
% [x, y] = parload( filename )
%
% input:
%--------
% filename: name for file
%
% output:
%--------
% x and y: variables contained in file
%
%==========================================================================


function [x, y] = parload(filename)


    foo = load(filename);
    whichVariables = fieldnames(foo);

    if (numel(whichVariables) == 1)
        x = foo.(whichVariables{1});
        y = [];
    elseif (numel(whichVariables) == 2)
        x = foo.(whichVariables{1});
        y = foo.(whichVariables{2});
    else
        error(['Problem with ', filename, '!']);
    end
    

end
