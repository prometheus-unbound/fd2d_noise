
function [parobj] = start_cluster(mode,folder_name,n_workers)
    
    
    if( strcmp(mode,'monch') )
        
        if( isempty(folder_name) )
            folder_name = getenv('SLURM_JOB_ID');
        end
        
        mkdir(folder_name);
        cluster = parallel.cluster.Generic('JobStorageLocation', folder_name);
        set(cluster, 'HasSharedFilesystem', true);
        set(cluster, 'ClusterMatlabRoot', '/apps/common/matlab/r2015a/');
        set(cluster, 'OperatingSystem', 'unix');
        set(cluster, 'IndependentSubmitFcn', @independentSubmitFcn);
        set(cluster, 'CommunicatingSubmitFcn', @communicatingSubmitFcn);
        set(cluster, 'GetJobStateFcn', @getJobStateFcn);
        set(cluster, 'DeleteJobFcn', @deleteJobFcn);
        
        parobj = parpool(cluster,n_workers,'IdleTimeout',120);
    
        
    elseif( strcmp(mode,'euler') || strcmp(mode,'brutus') )
        
        if( strcmp(mode,'euler') )
            cluster = parcluster('EulerLSF8h');
        elseif( strcmp(mode,'brutus') )
            cluster = parcluster('BrutusLSF8h');
        end
        
        if( isempty(folder_name) )
            folder_name = getenv('LSB_JOBID');
        end
        
        mkdir(folder_name);
        cluster.JobStorageLocation = folder_name;
        cluster.SubmitArguments = '-W 120:00 -R "rusage[mem=16384]"';
        parobj = parpool(cluster,16);
        
    end
    

end