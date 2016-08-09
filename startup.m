
fprintf( [ '\nIMPORTANT IMPORTANT:\n\n' ...
'When MATLAB tells you:\n' ...
'File ... is not found in the current folder or on the MATLAB path.\n\n' ...
'Please only click ''Change Folder'' and NOT on ''Add to Path''!\n\n\n' ] )

% cd tools
% addpath(genpath('../code'))
% addpath(genpath('../input'))
% addpath(genpath('../inversion'))
% addpath(genpath('../output'))
% addpath(genpath('../tools'))
% cd ..

current_path = pwd;
delimiter = current_path(1);
if( isempty( strfind(pwd, [delimiter 'fd2d_noise'] ) ) )
    error('\nPlease do not change the name of the parent folder ''fd2d_noise''.')
end


version_control = ver;

if( any( strcmpi( {version_control.Name}, 'matlab coder' ) ) )
    
    fprintf('You can use mex-functions!\n')
    fprintf('Try to compile wave propagation functions now...\n\n')
    cd code/mex_functions
    warning('off','all')
    compile_green
    compile_correlation
    % warning('on','all')
    cd ../..
    
    if( exit_code == 0 )
        fprintf('\nCompilation successful! Set use_mex to ''yes''!\n')
    else
        fprintf('\nCompilation NOT successful! Set use_mex to ''no''!\n')
    end
    
else
    
    fprintf('You do not have MATLAB Coder. Please set use_mex to ''no''.\n')
    
end

clear version_control
clear exit_code