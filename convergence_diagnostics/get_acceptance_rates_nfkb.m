%Code for Table S6, S7 NFkB with lasso

%We run 6 repeats for each of two starting parameter sets.
%If the 6 repeats have a PSRF < 1.2 for the energy chain then they are sampling the same energy basin and can be combined
%Here we combine the 6 repeats to get two groups and then we calculate
%average step acceptance rates and swap acceptance rates
addpath('../mcmc_diag');
%path to ptempest files and distributions
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
%path to results files
addpath('../SingleCellNFkB_pulse');
%%
mu = -25;
b = 2;
file_prefix = ['../SingleCellNFkB_pulse/SingleCellNFkB_continue_chain_reduced_model_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b)'];
c = 1;
LB = 5e4;
UB = 9.9e5;
repeat_index= [1:12];
for trajectory_number = [2,3,4]
    mu = -25;
    b = 2;
    num_repeats = length(repeat_index);
    num_pts1 = length(LB:UB);
    tmp = cell(1,num_repeats);
    dim = 26;
    par_matrix = zeros(num_pts1,dim,num_repeats);
    % Combine alternate repeats. Group 1 = [1,3,5,7,9,11]; Group 2 =
    % [2,4,6,8,10,12]
    for start = 1:2
        c = 1;
        sa = zeros(length(LB:UB),num_repeats/2);
        sw = zeros(length(LB:UB),num_repeats/2);
        for repeat = start:2:num_repeats
            display(['trajectory_number: ',num2str(trajectory_number),' repeat: ',num2str(repeat_index(repeat))]);
            data = load(['../SingleCellNFkB_pulse/SingleCellNFkB_continue_chain_reduced_model_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat_index(repeat)),'_progress990000.mat']);           
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
            sa(:,c) = data.step_acceptance(1,LB:UB);
            sw(:,c) = data.swap_acceptance(1,LB:UB);
            c=c+1;
        end
        mean_step_acceptance = sum(sum(sa))/((data.cfg.nsteps*length(LB:UB)*num_repeats/2))
        mean_swap_acceptance = sum(sum(sw))/((length(LB:UB)*num_repeats/2))
        %save(['nfkb_lasso_trajectory_',num2str(trajectory_number),'_group_',num2str(start),'_step_acceptance'],'mean_step_acceptance');
    end
end

