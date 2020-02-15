% path to parallel tempering scripts
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
% path to supplementary distributions
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
addpath('../lib')
addpath('../lib/aboxplot')
range = 5e4:1e3:9.9e5-1; %number of points to simulate to get likelihood distribution
e = zeros(2,length(range)*6,3);% array to store likelihood values. the first index is for with or without lasso. The second index corresponds to each parameter set tested. The third index corresponds to the NFkB trajectory
mu = -25;
b = 2;
repeat = 1;
counter=1;
trajectory_counter=1;
for trajectory_number = 2:4 % these are the trajectories for which the fits converged   
    c = 1; 
    energy_tmp = zeros(6,length(range)); % energy_tmp collects the likelihood values for a particular trajectory (all repeats)
    repeat_counter=1;
    %alternate repeats correspond to the same initial condition and are
    %combined
    for repeat = 1:2:12
        data = load(['../SingleCellNFkB_pulse/SingleCellNFkB_continue_chain_reduced_model_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat),'_progress990000.mat']);                     
        %Checking to make sure the initial parameter set is the same
        if repeat_counter==1
            tmp_init = data.cfg.initp;
        else
            if sum(sum(tmp_init~=data.cfg.initp))~=0
                disp('problem');
                return
            else 
                tmp_init = data.cfg.initp;
            end
        end   
        counter=1;
        for i = range
            i
            params = data.params_chain(1,:,i); 
            %calculate likelihood of each parameter set in the specified
            %range. in ptempest energy = liklihood - logprior. so now
            %likliehood  = energy + logprior
            energy_tmp(repeat_counter,counter) = data.energy_chain(1,i) + data.cfg.logpdf_prior_fcn(params);%likelihood(params,data1.cfg); 
            counter=counter+1;
        end
        repeat_counter=repeat_counter+1;
    end
    %Note that the (:) operatore concatenates columns. 
    %e(2,1,1) will correspond to energy_chain(1,range(1)) of trajectory 1, repeat 1
    %e(2,2,1) will correspond to energy_chain(1,range(1)) of trajectory 1, repeat 2
    %e(2,7,1) will correspond to energy_chain(1,range(2)) of trajectory 1, repeat 1
    e(2,:,trajectory_counter) = energy_tmp(:);
    %%
    energy_tmp1 = zeros(6,length(range));
    repeat_counter=1;
    for repeat = 1:2:12 % for biorxiv this was set to 2:2:12. Doesn't make a difference (since converged) but changed to 1:2:12 for consistency.
        repeat
        data2 = load(['../SingleCellNFkB_pulse/SingleCellNFkB_continue_chain_wolasso_fixedstart_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat),'_progress990000.mat']);       
        %Checking to make sure the initial parameter set is the same
        if repeat_counter==1
            tmp_init = data2.cfg.initp;
        else
            if sum(sum(tmp_init~=data2.cfg.initp))~=0
                disp('problem');
                return
            else 
                tmp_init = data2.cfg.initp;
            end
        end 
        % without lasso there is just a uniform prior and the energy = the
        % liklihood
        energy_tmp1(repeat_counter,:) = data2.energy_chain(1,range);
        repeat_counter=repeat_counter+1;
    end
    %Note that the (:) operatore concatenates columns. 
    %e(1,1,1) will correspond to energy_chain(1,range(1)) of trajectory 1, repeat 1
    %e(1,2,1) will correspond to energy_chain(1,range(1)) of trajectory 1, repeat 2.
    %e(1,7,1) will correspond to energy_chain(1,range(2)) of trajectory 1, repeat 1
    e(1,:,trajectory_counter) = energy_tmp1(:); 
    trajectory_counter = trajectory_counter+1;    
end
%%
f = figure;aboxplot(e,'OutlierMarker','none');
toPDF(f,'nfkb_updated_figures_Sept3/Figure4c_likelihoodcomparison_Sept2.pdf')
