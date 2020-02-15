%Code for Table S7
addpath('/shared2/LabUserFiles/Sanjana_Gupta/MoreNoise_of_model_reduction_project_automated/lib/mcmc_diag');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/model_reduction_project_automated/lib')
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
%%
mu = -25;
b = 2;
c = 1;
LB = 1e5;
UB = 5e5;
groups = {[1:8],[9:16]};
for trajectory_number =  1
    mu = -25;
    b = 2;
    num_repeats = 16;
    num_pts1 = length(LB:UB);
    tmp = cell(1,num_repeats);
    dim = 26;
    % Get first set of repeats
    for start = 1:2
        c = 1;
        sa = zeros(length(LB:UB),num_repeats/2);
        for repeat = groups{start}
            display(['trajectory_number: ',num2str(trajectory_number),' repeat: ',num2str(repeat)]);
            data = load(['/shared3/LabUserFiles/Sanjana_Gupta/LassoManuscript/SingleCellNFkB_reduced_model_continuous_adaptiveIKK_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat),'_progress500000.mat']);           
            if c==1
                tmp = data.cfg.initp;
            else
                if sum(sum(tmp~=data.cfg.initp))~=0
                    disp('problem');
                    return
                else 
                    tmp = data.cfg.initp;
                end
            end
            sa(:,c) = data.swap_acceptance(1,LB:UB);
            c=c+1;
        end
        mean_swap_acceptance = sum(sum(sa))/((length(LB:UB)*num_repeats/2))
        %save(['nfkb_lasso_trajectory_',num2str(trajectory_number),'_group_',num2str(start),'_swap_acceptance'],'mean_swap_acceptance');
    end
end

