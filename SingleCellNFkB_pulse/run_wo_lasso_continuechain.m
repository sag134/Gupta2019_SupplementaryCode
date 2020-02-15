% path to parallel tempering scripts
addpath('/shared2/LabUserFiles/Sanjana_Gupta/OriginalWithFixedStart/ptempest/core/');
% path to supplementary distributions
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
num_traj = 5;
b = 2;
%Location of previously run repeats that will be used to initialize the
%PT chains
file_prefix = 'initialize_wolasso/SingleCellNFkB_wolasso_randomstart_';

job = [1e29,1e29]; %dummy job index to initialize file
dlmwrite('jobs_wolassov2.txt',job,'-append');

for trajectory_number = 1:5
    for repeat = 1:12
        job = [trajectory_number,repeat];
        %open job file
        started_jobs = dlmread('jobs_wolassov2.txt');
        %check if this job has been started already
        [~,index] = ismember(job,started_jobs,'rows');
        if index==0
            dlmwrite('jobs_wolassov2.txt',job,'-append');
        else
            continue
        end
        if mod(repeat,2) == 1  
            full_prefix = [file_prefix,'trajectory_',num2str(trajectory_number),'_repeat_',num2str(1),'_progress'];
        else
            full_prefix = [file_prefix,'trajectory_',num2str(trajectory_number),'_repeat_',num2str(2),'_progress'];
        end
        S = dir([full_prefix,'*.*']);
        data = load(['initialize_wolasso/',S.name]);
        % load configuration file
        jobname = ['SingleCellNFkB_continue_chain_wolasso_fixedstart_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat)];
        cfg = config_wo_lasso_continue(jobname,trajectory_number);
        cfg.initp = data.params_chain(1,:,1e5);
        cfg.relstep_init= 0.001;
        cfg.max_init_steps = -1;
        % start parallel tempering
        parallel_tempering(cfg);  
    end
end
