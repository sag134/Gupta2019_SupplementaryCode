%We run 8 repeats for each of two starting parameter sets.
%If the 8 repeats have a PSRF < 1.2 for the energy chain then they are sampling the same energy basin and can be combined
%Here we combine the 6 repeats to get two groups and then we calculate
%parameter convergence using MPSRF/PSRF
addpath('/shared2/LabUserFiles/Sanjana_Gupta/MoreNoise_of_model_reduction_project_automated/mcmcdiag');
%path to ptempest files and distributions
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
%path to results files
%addpath('../SingleCellNFkB_pulse');
%%

groups = {[1:8],[9:16]};
mu = -25;
b = 2;
LB = 1e5;
UB =5e5;
max_num_pts = length(LB:UB);%Number of samples
num_repeats = 16;
for trajectory_number = [1]
    mu = -25;
    b = 2;
    tmp = cell(1,num_repeats);
    dim = 26; %Number of model parameters
    par_matrix = zeros(max_num_pts*num_repeats/2,dim,2);
    repeat_counter=1;
    for i = 1:2
        c=1;
        for repeat = groups{i}
            repeat
            display(['trajectory_number: ',num2str(trajectory_number),' repeat: ',num2str(repeat)]);
            data = load(['/shared3/LabUserFiles/Sanjana_Gupta/LassoManuscript/SingleCellNFkB_reduced_model_continuous_adaptiveIKK_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat),'_progress500000.mat']);           
            params = data.params_chain(1,:,LB:UB);
            num_pts = length(LB:UB);
            tmp{repeat_counter} = reshape(params,dim,num_pts); 
            par_matrix(c:c+num_pts-1,:,i) = tmp{repeat_counter}';
            c = c+num_pts;
            repeat_counter=repeat_counter+1;
        end
    end
    m = mpsrf(par_matrix)
    m1 = (psrf(par_matrix))
   %save(['nfkb_lasso_trajectory_',num2str(trajectory_number)],'m','m1','LB','UB');
end
