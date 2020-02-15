% path to parallel tempering scripts
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
% path to supplementary distributions
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
job = [1e29,1e29]; %dummy job index to initialize file
dlmwrite('jobs_wolasso.txt',job,'-append');
for trajectory_index = 1:5
    for repeat = 3:5       
        job = [trajectory_index,repeat];
        %open job file
        started_jobs = dlmread('jobs_wolasso.txt');
        %check if this job has been started already
        [~,index] = ismember(job,started_jobs,'rows');
        if index==0
            dlmwrite('jobs_wolasso.txt',job,'-append');
        else
            continue
        end
        % load configuration file
        jobname = ['SingleCellNFkB_wolasso_randomstart_trajectory_',num2str(trajectory_index),'_repeat_',num2str(repeat)];
        cfg = config_wo_lasso(jobname,trajectory_index);
        % start parallel tempering
        parallel_tempering(cfg);
    end
end

