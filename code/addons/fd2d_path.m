
%==========================================================================
% get path of fd2d_noise code
%
% [code_path] = fd2d_path()
%
%==========================================================================


function [code_path] = fd2d_path()


    current_path = pwd;
    folders = strsplit(current_path, filesep);
    id_project_folder = find(strncmp(folders, 'fd2d_noise', 10), 1, 'last');

    if (isempty(id_project_folder))
        error('Please stay within the ''fd2d_noise'' folder and keep the name ''fd2d_noise*'' (* = wildcard)!')
    end

    code_path = fullfile(filesep, folders{1:id_project_folder}, filesep);

    dir_content = ls(code_path);
    if (isempty(strfind(dir_content, 'code')))
        error('There is no code folder!')
    end
    if (isempty(strfind(dir_content, 'input')))
        error('There is no input folder!')
    end
    if (isempty(strfind(dir_content, 'inversion')))
        error('There is no inversion folder!')
    end
    if (isempty(strfind(dir_content, 'output')))
        error('There is no output folder!')
    end
    if (isempty(strfind(dir_content, 'tools')))
        error('There is no tools folder!')
    end


end
