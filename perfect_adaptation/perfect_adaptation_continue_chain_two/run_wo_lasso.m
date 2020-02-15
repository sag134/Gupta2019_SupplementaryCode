% this script will start-up a parallel tempering job for
%  the Michaelis-Menten model 
%
%   to execute from command line on remote server:
%   > nohup matlab -r run_pt > run_pt.out 2> run_pt.err < /dev/null &
  
% path to parallel tempering scripts
addpath('../../OriginalWithFixedStart/ptempest/core/');

% path to supplementary distributions
addpath('../../OriginalWithFixedStart/ptempest/core/distr/');
%%
% path to model-specific files (if not the current directory)
%addpath('~/googlecode/ptempest/examples/michment/');
for repeat = [1,2]
    data = load(['../perfect_adaptation/pa_wo_lasso_repeat_',num2str(repeat),'_progress900000.mat']);
    nchains = 4;
    jobname = ['perfect_adaptation_fixedstart_wo_lasso_repeat_',num2str(repeat)];
    %% load configuration file
    cfg = config_wo_lasso(nchains,jobname);
    cfg.initp = data.params_chain(1,:,end);
    cfg.max_init_steps = -1;
    cfg.relstep_init = 0.001;
    % start parallel tempering
    parallel_tempering(cfg);
end
