function communicatingSubmitFcn(cluster, job, props)
%COMMUNICATINGSUBMITFCN Submit a communicating MATLAB job to a SLURM cluster
%
% Set your cluster's CommunicatingSubmitFcn to this function using the following
% command:
%     set(cluster, 'CommunicatingSubmitFcn', @communicatingSubmitFcn);
%
% See also parallel.cluster.generic.communicatingDecodeFcn.
%

% Copyright 2010-2015 The MathWorks, Inc.

% Store the current filename for the errors, warnings and dctSchedulerMessages
currFilename = mfilename;
if ~isa(cluster, 'parallel.Cluster')
    error('parallelexamples:GenericSLURM:SubmitFcnError', ...
        'The function %s is for use with clusters created using the parcluster command.', currFilename)
end

decodeFunction = 'parallel.cluster.generic.communicatingDecodeFcn';

if ~cluster.HasSharedFilesystem
    error('parallelexamples:GenericSLURM:SubmitFcnError', ...
        'The submit function %s is for use with shared filesystems.', currFilename)
end

if ~strcmpi(cluster.OperatingSystem, 'unix')
    error('parallelexamples:GenericSLURM:SubmitFcnError', ...
        'The submit function %s only supports clusters with unix OS.', currFilename)
end

% The job specific environment variables
% Remove leading and trailing whitespace from the MATLAB arguments
matlabArguments = strtrim(props.MatlabArguments);
variables = {'MDCE_DECODE_FUNCTION', decodeFunction; ...
    'MDCE_STORAGE_CONSTRUCTOR', props.StorageConstructor; ...
    'MDCE_JOB_LOCATION', props.JobLocation; ...
    'MDCE_MATLAB_EXE', props.MatlabExecutable; ...
    'MDCE_MATLAB_ARGS', matlabArguments; ...
    'MDCE_DEBUG', 'true'; ...
    'MLM_WEB_LICENSE', props.UseMathworksHostedLicensing; ...
    'MLM_WEB_USER_CRED', props.UserToken; ...
    'MLM_WEB_ID', props.LicenseWebID; ...
    'MDCE_LICENSE_NUMBER', props.LicenseNumber; ...
    'MDCE_STORAGE_LOCATION', props.StorageLocation; ...
    'MDCE_CMR', cluster.ClusterMatlabRoot; ...
    'MDCE_TOTAL_TASKS', num2str(props.NumberOfTasks)};
% Set each environment variable to newValue if currentValue differs.
% We must do this particularly when newValue is an empty value,
% to be sure that we clear out old values from the environment.
for ii = 1:size(variables, 1)
    variableName = variables{ii,1};
    currentValue = getenv(variableName);
    newValue = variables{ii,2};
    if ~strcmp(currentValue, newValue)
        setenv(variableName, newValue);
    end
end

% Deduce the correct quote to use based on the OS of the current machine
if ispc
    quote = '"';
else
    quote = '''';
end


% The script name is communicatingJobWrapper.sh
scriptName = 'communicatingJobWrapper.sh';
% The wrapper script is in the same directory as this file
dirpart = fileparts(mfilename('fullpath'));
quotedScriptName = sprintf('%s%s%s', quote, fullfile(dirpart, scriptName), quote);

% Choose a file for the output. Please note that currently, JobStorageLocation refers
% to a directory on disk, but this may change in the future.
logFile = cluster.getLogLocation(job);
quotedLogFile = sprintf('%s%s%s', quote, logFile, quote);

jobName = sprintf('Job%d', job.ID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CUSTOMIZATION MAY BE REQUIRED %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% You might want to customize this section to match your cluster,
% for example to limit the number of nodes for a single job.
% You may also wish to supply additional submission arguments to
% the sbatch command here.
% % additionalSubmitArgs = sprintf('--partition=fichtner_compute_wk --time=07-00:00:00 --nodes=16 --ntasks-per-node=1 --mem-per-cpu=32768 --ntasks=%d', props.NumberOfTasks);
additionalSubmitArgs = sprintf('--partition=fichtner_compute_wk --time=07-00:00:00 --nodes=4 --ntasks-per-node=4 --mem-per-cpu=8192 --ntasks=%d', props.NumberOfTasks);
% additionalSubmitArgs = sprintf('--partition=fichtner_compute --time=00-05:00:00 --nodes=2 --ntasks-per-node=8 --mem-per-cpu=4096 --ntasks=%d', props.NumberOfTasks);
% additionalSubmitArgs = sprintf('--partition=fichtner_compute --time=01-00:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=32768 --ntasks=%d', props.NumberOfTasks);
% additionalSubmitArgs = sprintf('--partition=other_largemem --time=00-05:00:00 --nodes=2 --ntasks-per-node=8 --mem-per-cpu=8192 --ntasks=%d', props.NumberOfTasks);
% additionalSubmitArgs = sprintf('--partition=other_hugemem --time=01-00:00:00 --nodes=2 --ntasks-per-node=8 --mem-per-cpu=32768 --ntasks=%d', props.NumberOfTasks);
dctSchedulerMessage(5, '%s: Generating command for task %i', currFilename, ii);
commandToRun = getSubmitString(jobName, quotedLogFile, quotedScriptName, ...
    additionalSubmitArgs);

% Now ask the cluster to run the submission command
dctSchedulerMessage(4, '%s: Submitting job using command:\n\t%s', currFilename, commandToRun);
try
    % Make the shelled out call to run the command.
    [cmdFailed, cmdOut] = system(commandToRun);
catch err
    cmdFailed = true;
    cmdOut = err.message;
end
if cmdFailed
    error('parallelexamples:GenericSLURM:SubmissionFailed', ...
        'Submit failed with the following message:\n%s', cmdOut);
end

dctSchedulerMessage(1, '%s: Job output will be written to: %s\nSubmission output: %s\n', currFilename, logFile, cmdOut);

jobIDs = extractJobId(cmdOut);
% jobIDs must be a cell array
if isempty(jobIDs)
    warning('parallelexamples:GenericSLURM:FailedToParseSubmissionOutput', ...
        'Failed to parse the job identifier from the submission output: "%s"', ...
        cmdOut);
end
if ~iscell(jobIDs)
    jobIDs = {jobIDs};
end

% set the job ID on the job cluster data
cluster.setJobClusterData(job, struct('ClusterJobIDs', {jobIDs}));
