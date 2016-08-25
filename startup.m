

%- get current working directory ------------------------------------------
% current_path = pwd;


%- set path ---------------------------------------------------------------
addpath(genpath([pwd, filesep, 'code']))
addpath(genpath([pwd, filesep, 'input']))
addpath(genpath([pwd, filesep, 'inversion']))
addpath(genpath([pwd, filesep, 'output']))
addpath(genpath([pwd, filesep, 'tools']))


%- do not show warnings ---------------------------------------------------
warning off


%- for later runs, get path of fd2d_noise ---------------------------------
% folders = strsplit(current_path, filesep);
% id_project_folder = find(strncmp(folders, 'fd2d_noise', 10), 1, 'last');
% if (isempty(id_project_folder))
%     error('Please stay within the ''fd2d_noise'' folder and keep the name ''fd2d_noise*'' (* = wildcard)!')
% end
% fd2d_path = fullfile(filesep, folders{1:id_project_folder }, filesep);


%- check if mex functions can be used -------------------------------------
% version_control = ver;
% if( any( strcmpi( {version_control.Name}, 'matlab coder' ) ) )
%     fprintf('\nYou can use mex-functions!\n')
% end


%- cleanup ----------------------------------------------------------------
% clear current_path
% clear folders; clear id_project_folder
% clear fd2d_path; clear version_control
return



% if (any(strcmpi({version_control.Name}, 'matlab coder')))
%     
%     fprintf('Try to compile wave propagation functions now...\n\n')
%     cd([fd2d_path, 'code', filesep, 'mex_functions'])
%     
%     S = warning();
%     warning('off', 'all')
%     compile_green
%     compile_correlation
%     warning(S);
%     
%     cd(fd2d_path)
%     
%     if (exit_code == 0)
%         fprintf('\nCompilation successful! Set use_mex to ''yes''!\n\n')
%     else
%         fprintf('\nCompilation NOT successful! Set use_mex to ''no''!\n\n')
%     end
%     
% else
%     
%     fprintf('You do not have MATLAB Coder. Please set use_mex to ''no''.\n\n')
%     
% end
% 
% %- cleanup ----------------------------------------------------------------
% clear folders; clear id_project_folder
% clear fd2d_path; clear version_control
% clear S


