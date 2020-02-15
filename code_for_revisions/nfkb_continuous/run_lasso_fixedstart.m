% path to parallel tempering scripts
addpath('/shared2/LabUserFiles/Sanjana_Gupta/OriginalWithFixedStart/ptempest/core/');
% path to supplementary distributions
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');

B_VALUES = [2];
job = [1e29,1e29,1e29,1e29]; %dummy job index to initialize file
dlmwrite('jobs_lasso1.txt',job,'-append');
%initializing from a previous run where b was very large and therefore easy
%to fit
data = load('SingleCellNFkB_reduced_model_continuous_adaptiveIKK_randomstart_withlasso_mu_-25_b_10_trajectory_1_repeat_1_progress300000.mat');
for mu = [-25]
    for b = B_VALUES
        for trajectory_number = 1
            for repeat = 9:12%[1:8]%1:12
                job = [mu,b,trajectory_number,repeat];
                %open job file
                started_jobs = dlmread('jobs_lasso1.txt');
                %check if this job has been started already
                [~,index] = ismember(job,started_jobs,'rows');
                if index==0
                    dlmwrite('jobs_lasso1.txt',job,'-append');
                else
                    continue
                end
                % load configuration file
                jobname = ['SingleCellNFkB_reduced_model_continuous_adaptiveIKK_randomstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat)];
                cfg = config_lasso(jobname,mu,b,trajectory_number);
                %cfg.initp = data.params_chain(1,:,3e5);
                if repeat<11
                    cfg.initp = data.params_chain(1,:,2.8e5);
                else
                    cfg.initp = data.params_chain(1,:,3e5);
                end
                cfg.relstep_init = 0.001;
                cfg.max_init_steps = -1;
                % start parallel tempering
                parallel_tempering(cfg);
            end
        end
    end
end
