% path to parallel tempering scripts
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
% path to supplementary distributions
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
%% 
nchains = 4;
MU_VALUES = [-10]%,-8];
B_VALUES = 1%[0.1,0.2,0.5,1,5,0.01,0.05];
job = [1e29,1e29,1e29,1e29]; %dummy job index to initialize file
dlmwrite('jobs_lasso.txt',job,'-append');
for fraction_error=[0.05]%,0.1]
    for mu = MU_VALUES
        for b = B_VALUES
            for repeat = 1:4
                job = [mu,b,fraction_error,repeat];
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
                jobname = ['five_node_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_repeat_',num2str(repeat)];
                cfg = config_lasso(nchains,jobname,mu,b,fraction_error);
                % start parallel tempering
                parallel_tempering(cfg);
            end
        end
    end
end