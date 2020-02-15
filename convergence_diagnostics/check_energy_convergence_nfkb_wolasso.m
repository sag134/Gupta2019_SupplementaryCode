%Code for Table S4 without lasso

%We run 6 repeats for each of two starting parameter sets.
%If the 6 repeats have a PSRF < 1.2 for the energy chain then they are sampling the same energy basin and can be combined
%Third party code used for psrf / mpsrf calculations available here: https://research.cs.aalto.fi/pml/software/mcmcdiag/
addpath('../mcmc_diag');
%path to ptempest files and distributions
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
%path to results files
addpath('../SingleCellNFkB_pulse');
%%
c = 1;
%range of chain to consider
LB = 5e4;
UB = 9.9e5;
repeat_index= [1:12];
for trajectory_number = [2,3,4] %For the paper we used trajectories 2,3,4 as 1,5 had some convergence issues even though the regularization result was qualitatively the same
    mu = -25;
    b = 2;
    num_repeats = length(repeat_index);
    num_pts1 = length(LB:UB);
    tmp = cell(1,num_repeats);
    dim = 26;
    par_matrix = zeros(num_pts1,dim,num_repeats);
    %Alternate repeats e.g. [1,3,5 ..] come from the same initial parameter
    %set and are checked for energy convergence
    for start = 1:2
        figure
        c = 1;
        e = zeros(length(LB:UB),1,num_repeats/2); %The shape of e matches the input requirement for the psrf function. an NxDxM matrix where N is the number of samples, D is the dimension, and M is the number of MCMC simulations. here the dimension is 1.
        for repeat = start:2:num_repeats
            subplot(2,3,c);
            display(['trajectory_number: ',num2str(trajectory_number),' repeat: ',num2str(repeat_index(repeat))]);
            data = load(['../SingleCellNFkB_pulse/SingleCellNFkB_continue_chain_wolasso_fixedstart_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat_index(repeat)),'_progress990000.mat']);           
            %Checking to make sure the initial parameter set is the same
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
            %Collect the energy chain
            e(:,1,c) = data.energy_chain(1,LB:UB);
            histogram(e(:,1,c),100,'EdgeColor','none');
            c=c+1;
        end
        energy_psrf = psrf(e);
        save(['nfkb_wolasso_trajectory_',num2str(trajectory_number),'_group_',num2str(start),'_energy_psrf'],'energy_psrf','LB','UB');
    end
end
