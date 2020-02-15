%Code for Table S6
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
    num_repeats = 16
    num_pts1 = length(LB:UB);
    tmp = cell(1,num_repeats);
    dim = 26;
    par_matrix = zeros(num_pts1,dim,num_repeats);
    % Get first set of repeats
    for start = 1:2
        c = 1;
        sa = zeros(length(LB:UB),num_repeats/2);
        sw = zeros(length(LB:UB),num_repeats/2);
        for repeat = groups{start}
            display(['trajectory_number: ',num2str(trajectory_number),' repeat: ',num2str(repeat)]);
            data = load(['/shared3/LabUserFiles/Sanjana_Gupta/LassoManuscript/SingleCellNFkB_reduced_model_continuous_adaptiveIKK_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat),'_progress500000.mat']);           
            if c==1
                tmp = data.cfg.initp;
            else
                if sum(sum(tmp~=data.cfg.initp))~=0 %Make sure that repeats in the same group have the same initial conditions.
                    disp('problem');
                    return
                else 
                    tmp = data.cfg.initp;
                end
            end
            sa(:,c) = data.step_acceptance(1,LB:UB);
            sw(:,c) = data.swap_acceptance(1,LB:UB);
            c=c+1;
        end
        mean_step_acceptance = sum(sum(sa))/((data.cfg.nsteps*length(LB:UB)*num_repeats/2))
        mean_swap_acceptance = sum(sum(sw))/((length(LB:UB)*num_repeats/2))
        %save(['nfkb_lasso_trajectory_',num2str(trajectory_number),'_group_',num2str(start),'_step_acceptance'],'mean_step_acceptance','mean_swap_acceptance','sa','sw');
    end
end

