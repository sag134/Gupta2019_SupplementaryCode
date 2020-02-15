% path to parallel tempering scripts. OriginalWithFixedStart contains a
% modified chain initialization file to enable intializing the runs to the
% parameter set in variable cfg.initp

addpath('/shared2/LabUserFiles/Sanjana_Gupta/OriginalWithFixedStart/ptempest/core/');
% path to supplementary distributions
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');

%A shorter chain fit to one of the trajectories is used to initialize fits
%to all the trajectories
data = load('SingleCellNFkB_reduced_model_randomstart_withlasso_mu_-22_b_2_trajectory_1_repeat_1_progress290000.mat');

B_VALUES = 2;
job = [1e29,1e29,1e29,1e29]; %dummy job index to initialize file
dlmwrite('jobs_lasso.txt',job,'-append');
for mu = -25
    for b = B_VALUES
        for trajectory_number = 1:5
            for repeat = 1:2
                job = [mu,b,trajectory_number,repeat];
                %open job file
                started_jobs = dlmread('jobs_lasso.txt');
                %check if this job has been started already
                [~,index] = ismember(job,started_jobs,'rows');
                if index==0
                    dlmwrite('jobs_lasso.txt',job,'-append');
                else
                    continue
                end
                % load configuration file
                jobname = ['SingleCellNFkB_reduced_model_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat)];  
                cfg = config_lasso(jobname,mu,b,trajectory_number);
                cfg.initp = data.params_chain(1,:,2e5);
                cfg.relstep_init= 0.001;
                cfg.max_init_steps = -1;
                % start parallel tempering
                parallel_tempering(cfg);
            end
        end
    end
end
