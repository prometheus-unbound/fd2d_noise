

%- do not show warnings ---------------------------------------------------
warning off


%- set path ---------------------------------------------------------------
addpath(genpath([pwd, filesep, 'code']))
addpath(genpath([pwd, filesep, 'input']))
addpath(genpath([pwd, filesep, 'inversion']))
addpath(genpath([pwd, filesep, 'literature']))
addpath(genpath([pwd, filesep, 'output']))
addpath(genpath([pwd, filesep, 'tools']))


%- print message ----------------------------------------------------------
fprintf(['\nPath is set correctly!' ... 
        '\nYou can start with the tutorial!\n\n']);