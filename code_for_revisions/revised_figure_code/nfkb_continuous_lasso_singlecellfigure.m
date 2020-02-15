addpath('../lib')
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
addpath('/shared3/LabUserFiles/Sanjana_Gupta/LassoManuscript/');
%%
for group_number=1:2
    step = 2e4; % step interval at which to simulate single cell trajectories
    f = figure('units','normalized','outerposition',[0,0,1,0.5]);
    c = 1;
    if group_number==1
        repeat_index =1:8;
    else
        repeat_index =9:16;
    end
    for trajectory_number = 1
        mu = -25;
        b = 2;
        num_repeats = length(repeat_index);   
        tmp = cell(1,num_repeats);
        dim = 26;
        LB = 1e5;
        UB = 5e5;
        num_pts1 = length(LB:UB); % num_pts1 * 6 is the total number of points in each histogram
        par_matrix = zeros(num_pts1,dim,num_repeats);
        for repeat = 1:num_repeats
            display(['trajectory_number: ',num2str(trajectory_number),' repeat: ',num2str(repeat)]);
            %alternate repeats correspond to the same initial condition and are
            %combined
            data = load(['../../LassoManuscript/SingleCellNFkB_reduced_model_continuous_adaptiveIKK_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat_index(repeat)),'_progress500000.mat']);           
            %Checking to make sure the initial parameter set is the same
            %within group but not between groups
            if repeat_index(repeat) == 9
              if sum(sum(tmp_init~=data.cfg.initp))==0
                    disp('problem');
                    return
              end  
            end
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
            %collect parameters from this repeat
            par_matrix(:,:,repeat) = tmp{repeat}';
            %plot single-cell fits
            for j = LB:step:UB % length(LB:step:UB) * num_repeats is the total number of trajectories plotted = 288 here
                subplot(3,4,c);
                params = data.params_chain(1,:,j);
                [err, sp, obsv] = simulate( [],[], params,data.cfg);
                p = plot(obsv(1:end-1),'b','LineWidth',2);
                p.Color(4) = 0.1;
                hold on
                ylim([0.5 4])
            end
        end
        hold on
        errorbar(1:length(data.cfg.data{1}.mean)-1,data.cfg.data{1}.mean(1:end-1),data.cfg.data{1}.stdev(1:end-1),'r','LineWidth',1)
        c=c+1;
        %% plot the group lasso regularization parameters
        tmp = 0;
        for j = [24,25,26]
            cfg=data.cfg;
            subplot(3,4,c);c=c+1; 
            params = par_matrix(:,j,:); %extract the relevent parameter value
            h = histogram(params(:),1000,'EdgeColor','none');
            mxhst = max(h.Values);
            hold on
            ylim([0 mxhst]) %scale y axis to max value of distribution
            % plot the laplace equation scaled so that the max value of the
            % laplace is the max value of the distribution (since we are
            % plotting frequency distributions and not probability
            % distributions so we don't have to integrate to 1)
            xlim([cfg.param_defs{j}.min,cfg.param_defs{j}.max])
            l = @(x,b,mu)(1/(2*b))*exp(-abs(x-mu)/b); %laplace equation 
            yl = ylim; %yl(2) is the desired max value of the laplace
            xlim([cfg.param_defs{j}.min,cfg.param_defs{j}.max])
            mx  = (1/(2*b)); %max value of laplace without scaling
            xl = xlim;
            ll = yl(2)*l(xl(1):0.0005:xl(2),b,mu)/(mx);
            plot(xl(1):0.0005:xl(2),ll,'m','LineWidth',2);
            set(gca,'YTickLabels',[]);
        end
        %return
       % fig_print(f,['GroupLasso_trajectory_',num2str(trajectory_number),'.png'],[5000,1000],1,1,'w');
    end
   % toPDF(f,['../RevisedFiguresMATLAB/nfkb_continuous_group_lasso_fit_figure_group_',num2str(group_number),'Dec13.pdf']);
    return
end





    