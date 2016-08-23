
function check_mex_files(use_mex)


    if (strcmp(use_mex, 'no'))

        if (exist([fd2d_path(), 'code', filesep, 'mex_functions', filesep, 'codegen'], 'dir'))
            S = warning();
            warning('off', 'all')
            rmdir([fd2d_path(), 'code', filesep, 'mex_functions', filesep, 'codegen'], 's')
            warning(S);
        end
        delete([fd2d_path(), 'code', filesep, 'mex_functions', filesep, 'run*'])
        copyfile([fd2d_path(), 'code', filesep, 'run1_forward_green.m'], [fd2d_path(), 'code', filesep, 'mex_functions', filesep, 'run1_forward_green_mex.m'])
        copyfile([fd2d_path(), 'code', filesep, 'run2_forward_correlation.m'], [fd2d_path(), 'code', filesep, 'mex_functions', filesep, 'run2_forward_correlation_mex.m'])
        copyfile([fd2d_path(), 'code', filesep, 'run3_adjoint.m'], [fd2d_path(), 'code', filesep, 'mex_functions', filesep, 'run3_adjoint_mex.m'])

    else

        delete([fd2d_path(), 'code', filesep, 'mex_functions', filesep, 'run*.m'])

        dir_content = ls([fd2d_path(), 'code', filesep, 'mex_functions']);
        if (isempty(strfind(dir_content, 'run1_forward_green')))
            warning('run1_forward_green.m was not compiled! Try running startup.m again!')
            copyfile([fd2d_path(), 'code', filesep, 'run1_forward_green.m'], [fd2d_path(), 'code', filesep, 'mex_functions', filesep, 'run1_forward_green_mex.m'])
        end

        if (isempty(strfind(dir_content, 'run2_forward_correlation')))
            warning('run2_forward_correlation.m was not compiled! Try running startup.m again!')
            copyfile([fd2d_path(), 'code', filesep, 'run2_forward_correlation.m'], [fd2d_path(), 'code', filesep, 'mex_functions', filesep, 'run2_forward_correlation_mex.m'])
        end

        if (isempty(strfind(dir_content, 'run3_adjoint')))
            warning('run3_adjoint.m was not compiled! Try running startup.m again!')
            copyfile([fd2d_path(), 'code', filesep, 'run3_adjoint.m'], [fd2d_path(), 'code', filesep, 'mex_functions', filesep, 'run3_adjoint_mex.m'])
        end

    end


end
