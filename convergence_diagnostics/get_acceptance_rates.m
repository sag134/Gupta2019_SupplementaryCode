%Code for synthetic examples in Table S6, S7
%------------------------------------------
%Some details about location of results, hyperparameters to be used, and
%start and end points of the chains to be considered. "Experiment" refers
%to each example
experiments = { ...
    struct('name','three_node','mu',-10,'b',1,'LB',5e5,'UB',9e5,'num_par',6) ...
    struct('name','five_node','mu',-10,'b',1,'LB',5e5,'UB',9e5,'num_par',20)...
    struct('name','linear_full_trajectory','mu',-10,'b',0.5,'LB',7e5,'UB',9e5,'num_par',6)...
    struct('name','perfect_adaptation_fixedstart','mu',-10,'b',1,'LB',1e5,'UB',9e5,'num_par',6,'subfolder','perfect_adaptation_continue_chain_two')...
    };
nrepeats = 2; %report results for two repeats
nexpt = length(experiments);
step_acceptance_rate_lasso = zeros(nexpt,nrepeats);
step_acceptance_rate_wo_lasso = zeros(nexpt,nrepeats);
swap_acceptance_rate_lasso = zeros(nexpt,nrepeats);
swap_acceptance_rate_wo_lasso = zeros(nexpt,nrepeats);
for i = 1:nexpt
    folder_name = experiments{i}.name;
    if isfield(experiments{i},'subfolder') %check if results are in a subfolder or in top level directory
        folder_name = [folder_name,'/',experiments{i}.subfolder];
    end
    expt_name = experiments{i}.name;
    mu = experiments{i}.mu;
    b = experiments{i}.b;
    LB = experiments{i}.LB;
    UB = experiments{i}.UB;
    for repeat = 1:nrepeats
        % WITH LASSO
        data1 = load(['../',folder_name,'/',expt_name,'_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_repeat_',num2str(repeat),'_progress',num2str(UB),'.mat']);
        % calculate average step acceptance rate: number of accepted steps
        % / total number of steps; total number of steps = number of swaps
        % * number of steps per swap
        step_acceptance_rate_lasso(i,repeat) = sum(data1.step_acceptance(1,LB:UB))/(data1.cfg.nsteps*(UB-LB+1));
        % calculate average swap acceptance rate: number of accepted swaps
        % / total number of swaps; 
        swap_acceptance_rate_lasso(i,repeat) = sum(data1.swap_acceptance(1,LB:UB))/((UB-LB+1));
        % WITHOUT LASSO
        data2 = load(['../',folder_name,'/',expt_name,'_wo_lasso_repeat_',num2str(repeat),'_progress',num2str(UB),'.mat']);
        step_acceptance_rate_wo_lasso(i,repeat) = sum(data2.step_acceptance(1,LB:UB))/(data2.cfg.nsteps*(UB-LB+1));
        swap_acceptance_rate_wo_lasso(i,repeat) = sum(data2.swap_acceptance(1,LB:UB))/((UB-LB+1));
    end
end
