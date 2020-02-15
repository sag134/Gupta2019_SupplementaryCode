% path to parallel tempering scripts
addpath('../../OriginalWithFixedStart/ptempest/core/');
% path to supplementary distributions
addpath('../../OriginalWithFixedStart/ptempest/core/distr/');
%% 

nchains = 4;
MU_VALUES = [-10,-8];
B_VALUES = [5,1,0.5,0.1];%,0.05,0.01];
job = [1e29,1e29,1e29]; %dummy job index to initialize file
dlmwrite('jobs_lasso.txt',job,'-append');
for mu = sort(MU_VALUES)
    for b = B_VALUES        
        for repeat = 1:2
            data = load(['../perfect_adaptation_continue_chain/perfect_adaptation_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_repeat_',num2str(repeat),'_progress900000.mat']);
            job = [mu,b,repeat];
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
            jobname = ['perfect_adaptation_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_repeat_',num2str(repeat)];
            cfg = config_lasso(nchains,jobname,mu,b);            
            cfg.initp = data.params_chain(1,:,end);
            cfg.max_init_steps = -1;
            cfg.relstep_init = 0.001;
            % start parallel tempering
            parallel_tempering(cfg);
            quit;
        end
    end
end