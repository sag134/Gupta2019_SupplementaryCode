%Code for Table S2,S3 NFkB with lasso

%We run 6 repeats for each of two starting parameter sets.
%If the 6 repeats have a PSRF < 1.2 for the energy chain then they are sampling the same energy basin and can be combined
%Here we combine the 6 repeats to get two groups and then we calculate
%parameter convergence using MPSRF/PSRF
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
LB = 5e4;
UB =9.9e5;
max_num_pts = length(LB:UB);%Number of samples
repeat_index= 1:12;
for trajectory_number = [2,3,4]
    mu = -25;
    b = 2;
    num_repeats = length(repeat_index);
    tmp = cell(1,num_repeats);
    dim = 26; %Number of model parameters
    par_matrix = zeros(max_num_pts*num_repeats/2,dim,2);
    % Combine alternate repeats. Group 1 = [1,3,5,7,9,11]; Group 2 =
    % [2,4,6,8,10,12]
    for i = 1:2
        c=1;
        for repeat = i:2:num_repeats
            display(['trajectory_number: ',num2str(trajectory_number),' repeat: ',num2str(repeat_index(repeat))]);
            data = load(['../SingleCellNFkB_pulse/SingleCellNFkB_continue_chain_reduced_model_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat_index(repeat)),'_progress990000.mat']);           
            params = data.params_chain(1,:,LB:UB);
            num_pts = length(LB:UB);
            tmp{repeat} = reshape(params,dim,num_pts); 
            par_matrix(c:c+num_pts-1,:,i) = tmp{repeat}';
            c = c+num_pts;
        end
    end
    m = mpsrf(par_matrix)
    m1 = max(psrf(par_matrix))
%   save(['nfkb_lasso_trajectory_',num2str(trajectory_number)],'m','m1','LB','UB');
end
