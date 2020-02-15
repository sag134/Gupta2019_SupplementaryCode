% path to parallel tempering scripts
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
% path to supplementary distributions
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');

B_VALUES = [2];%[0.1 0.2,0.5,1,2];
job = [1e29,1e29,1e29,1e29]; %dummy job index to initialize file
dlmwrite('jobs_lasso.txt',job,'-append');
for mu = [-25]%,-20]
    for b = B_VALUES
        for trajectory_number = 1%:5
            for repeat = 1:10
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
                jobname = ['SingleCellNFkB_reduced_model_continuous_adaptiveIKK_randomstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat)];
                cfg = config_lasso(jobname,mu,b,trajectory_number);
                % start parallel tempering
                parallel_tempering(cfg);
            end
        end
    end
end
