%Code for synthetic examples in Tables S2 and S3
%-----------------------------------------------
%Third party code used for psrf / mpsrf calculations available here: https://research.cs.aalto.fi/pml/software/mcmcdiag/
addpath('../mcmc_diag');
%path to ptempest files and distributions
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');

folder_path = './convergence_test_results/';

%Some details about location of results, hyperparameters to be used, and
%start and end points of the chains to be considered
experiments = { ...
    struct('name','three_node','mu',-10,'b',1,'LB',5e5,'UB',9e5,'num_par',6) ...
    struct('name','five_node','mu',-10,'b',1,'LB',5e5,'UB',9e5,'num_par',20)...
    struct('name','linear_full_trajectory','mu',-10,'b',0.5,'LB',7e5,'UB',9e5,'num_par',6)...
    struct('name','perfect_adaptation_fixedstart','mu',-10,'b',1,'LB',1e5,'UB',9e5,'num_par',6,'subfolder','perfect_adaptation_continue_chain_two')...
    };

nrepeats = 2; % consider 2 repeats per example

nexpt = length(experiments);
lasso_mpsrf = zeros(1,nexpt); 
wo_lasso_mpsrf = zeros(1,nexpt);
lasso_psrf = zeros(1,nexpt); 
wo_lasso_psrf = zeros(1,nexpt);

for i =  1:nexpt
    folder_name = experiments{i}.name;
    if isfield(experiments{i},'subfolder') %check if results are in a subfolder or in top level directory
        folder_name = [folder_name,'/',experiments{i}.subfolder];
    end
    
    addpath(['../',folder_name])
    expt_name = experiments{i}.name;
    mu = experiments{i}.mu;
    b = experiments{i}.b;
    LB = experiments{i}.LB;
    UB = experiments{i}.UB;
    num_par = experiments{i}.num_par;
    num_pts = UB-LB+1;
    tmp = cell(1,nrepeats);
    par_matrix = zeros(num_pts,num_par,nrepeats);
    
    % WITH LASSO
    for repeat = 1:nrepeats
        %load results file
        data = load(['../',folder_name,'/',expt_name,'_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_repeat_',num2str(repeat),'_progress',num2str(UB),'.mat']);
        % get parameter chain
        params = data.params_chain(1,1:num_par,LB:UB);
        % save parameter chain from each repeat
        tmp{repeat} = reshape(params,num_par,num_pts);
        par_matrix(:,:,repeat) = tmp{repeat}';
    end
    %calculate mpsrf
    lasso_mpsrf(i) = mpsrf(par_matrix);
    %calculate the max psrf value across all parameters
    lasso_psrf(i) = max(psrf(par_matrix));
    %WITHOUT LASSO
    for repeat = 1:nrepeats
        %load results file
        data = load(['../',folder_name,'/',expt_name,'_wo_lasso_repeat_',num2str(repeat),'_progress',num2str(UB),'.mat']);
        % get parameter chain
        params = data.params_chain(1,1:num_par,LB:UB);
        tmp{repeat} = reshape(params,num_par,num_pts);
        % save parameter chain from each repeat
        par_matrix(:,:,repeat) = tmp{repeat}';
    end
    %calculate mpsrf
    wo_lasso_mpsrf(i) = mpsrf(par_matrix);
    %calculate the max psrf value across all parameters
    wo_lasso_psrf(i) =  max(psrf(par_matrix));
end

%save([folder_path,example,'_wo_lasso_psrf'],'par_matrix','tmp','wo_lasso_psrf')

