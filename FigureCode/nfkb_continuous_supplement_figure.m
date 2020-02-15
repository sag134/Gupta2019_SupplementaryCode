rmpath('../SingleCellNFkB_pulse')
addpath('/shared2/LabUserFiles/Sanjana_Gupta/MoreNoise_of_model_reduction_project_automated/lib/mcmc_diag');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/model_reduction_project_automated/lib')
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
addpath('../SingleCellNFkB_continuous');
addpath('../lib/');

step = 2e4;
repeat_index = 2:5;
for trajectory_number = 1
    f = figure('units','normalized','outerposition',[0,0,1,0.5]);
    c = 1;
    mu = -25;
    b = 2;
    num_repeats = length(repeat_index);   
    tmp = cell(1,num_repeats);
    dim = 26;
    LB = 5e4;
    UB = 9e5;
    num_pts1 = length(LB:UB);
    par_matrix = zeros(num_pts1,dim,num_repeats);
    % Get first set of repeats
    for repeat = 1:num_repeats
        display(['trajectory_number: ',num2str(trajectory_number),' repeat: ',num2str(repeat_index(repeat))]);
        data = load(['../SingleCellNFkB_continuous/SingleCellNFkB_reduced_model_continuous_adaptiveIKK_randomstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat_index(repeat)),'_progress900000.mat']);           
        
        %Checking to make sure the initial parameter set is the same
        if repeat==1
            tmp_init = data.cfg.initp;
        else
            if sum(sum(tmp_init~=data.cfg.initp))~=0
                disp('problem');
                return
            else 
                tmp_init = data.cfg.initp;
            end
        end   
        
        params = data.params_chain(1,:,LB:UB);
        tmp{repeat} = reshape(params,dim,num_pts1);
        par_matrix(:,:,repeat) = tmp{repeat}';
        for j = LB:step:UB
            subplot(1,4,c);
            params = data.params_chain(1,:,j);
            [err, sp, obsv] = simulate( [],[], params,data.cfg);
            if err==1
                display('error')
                return
            end
            p = plot(obsv(1:end-1),'b','LineWidth',2);
            p.Color(4) = 0.1;
            hold on
            ylim([0.5 5])
        end
    end
    hold on
    errorbar(1:length(data.cfg.data{1}.mean)-1,data.cfg.data{1}.mean(1:end-1),data.cfg.data{1}.stdev(1:end-1),'r','LineWidth',1)
    c=c+1;
    %%
    tmp = 0;
    for j = [24,25,26]
        cfg=data.cfg;
        subplot(1,4,c);c=c+1; 
        params = par_matrix(:,j,:);
        h = histogram(params(:),1000,'EdgeColor','none','FaceAlpha',0.7);
        mxhst = max(h.Values);
        hold on
        if mxhst>tmp
            tmp = mxhst;
            ylim([0 mxhst])
        end
        yl = ylim;
        plot([mu,mu],yl,'--r','LineWidth',2);
        xlim([cfg.param_defs{j}.min,cfg.param_defs{j}.max])

        l = @(x,b,mu)(1/(2*b))*exp(-abs(x-mu)/b); %laplace equation 
        yl = ylim;
        xlim([cfg.param_defs{j}.min,cfg.param_defs{j}.max])
        mx  = (1/(2*b));
        xl = xlim;
        ll = yl(2)*l(xl(1):0.0005:xl(2),b,mu)/(mx);
        plot(xl(1):0.0005:xl(2),ll,'m','LineWidth',2);
        set(gca,'YTickLabels',[]);
    end
    %%
    %fig_print(f,['GroupLasso_continuous_treatment.pdf'],[5000,1000],1,1,'w');
    %savefig(['GroupLasso_continuous_treatment.pdf'],'pdf',f);
    %toPDF(f,'GroupLasso_continuous_treatment.pdf')
end







    