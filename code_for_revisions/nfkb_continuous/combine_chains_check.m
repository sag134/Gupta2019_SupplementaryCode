%Code for Table S5, checking energy convergence before combining independent chains
addpath('/shared2/LabUserFiles/Sanjana_Gupta/MoreNoise_of_model_reduction_project_automated/lib/mcmc_diag');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/model_reduction_project_automated/lib')
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
addpath('../SingleCellNFkB_pulse');
%%
mu = -25;
b = 2;
file_prefix = ['SingleCellNFkB_reduced_model_continuous_adaptiveIKK_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b)];
c = 1;
groups = {[1:8],[9:16]};
LB = 1e5;
UB = 5e5;
for trajectory_number = 1
    mu = -25;
    b = 2;
    num_repeats = 8;
    num_pts1 = length(LB:UB);
    tmp = cell(1,num_repeats);
    dim = 26;
    par_matrix = zeros(num_pts1,dim,num_repeats);
    % Get first set of repeats
    for start = 1:2
        figure
        c = 1;
        e = zeros(length(LB:UB),1,num_repeats);
        for repeat = groups{start}
            display(['trajectory_number: ',num2str(trajectory_number),' repeat: ',num2str(repeat)]);
            data = load(['/shared3/LabUserFiles/Sanjana_Gupta/LassoManuscript/SingleCellNFkB_reduced_model_continuous_adaptiveIKK_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat),'_progress500000.mat']);           
            if c==1
                tmp = data.cfg.initp;
            else
                if sum(sum(tmp~=data.cfg.initp))~=0 %Check to make sure that the start points within a group are the same
                    disp('problem');
                    return
                else 
                    tmp = data.cfg.initp;
                end
            end
            e(:,1,c) = data.energy_chain(1,LB:UB);
            c=c+1;
        end
        energy_psrf = psrf(e)
       %save(['nfkb_continuous_lasso_trajectory_',num2str(trajectory_number),'_group_',num2str(start),'_energy_psrf'],'energy_psrf');
    end
end

