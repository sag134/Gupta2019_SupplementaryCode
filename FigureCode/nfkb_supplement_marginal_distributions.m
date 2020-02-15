addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
addpath('../lib');
%%
fs = 12;
%Group lasso penalty index for each rate constant parameter
multiplier = [25,25,25,25,26,26,25,25,26,26,26,26,26,25,25,26,24,24,24,24];
%parameters (proteins and rate constants) from Pekalski et al. 2013 (model units were all in terms of "#
%molecules". no conversions necessary).
L_par = log10([1e5,2e5,7000,1.2e-5,1.2e-3,nan,nan,5e-7,nan,1e-2,nan,2e-3,5e-3,5e-2,1e-4,2e-5,nan,nan]);
%parameters (proteins and rate constants)  from Lee et al. 2014
conv = 0.04/50000; % conversion factor from micromolar to molecules units from fold change paper
R_par = log10([50000,nan,nan,nan,nan,nan,nan,0.5*conv,0.05,0.0026,0.000052,0.00067,0.000335,0.01,0.0005,0.000022,nan,nan,nan]);
%parameters from Kearns et al. 2006. Model units are in micromolar so same
%conversion factor is used. Also minutes is converted into seconds
W_par = log10([0.0875/conv,0.8/conv,nan,nan,nan,nan,nan,30 * conv / 60, 0.00006 / 60,5.4/60,0.0048/60,0.018/60,0.012/60,0.828/60,0.12/60,0.00006/60,nan,nan,nan]);
clr_fr = [0.5,1,1];
tcounter=0;
repeat_index = 1:2:12;
for trajectory_number = 2 % showing results for a representative trajectory
    tcounter=tcounter+1;
    LB = 5e4;
    UB = 9.9e5;
    num_repeats = length(repeat_index);
    dim = 26;
    num_pts = UB-LB+1;
    par_matrix = zeros(num_pts,dim,num_repeats);
    tmp = cell(1,num_repeats);
    % get data
    mu = -25;
    b = 2;
    emin = 1e29;
    par_min = zeros(1,dim);
    for repeat = 1:num_repeats
        display(['trajectory_number: ',num2str(trajectory_number),' repeat: ',num2str(repeat_index(repeat))]);
        data = load(['../SingleCellNFkB_pulse/SingleCellNFkB_continue_chain_reduced_model_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat_index(repeat)),'_progress990000.mat']);    
        e = data.energy_chain(1,:);  
        % From the first repeat get a representative best-fit parameter set
        if repeat==1
            ind = find(e ==min(e),1);
            emin = e(ind);
            par_min = data.params_chain(1,:,ind);
        end     
        params = data.params_chain(1,:,LB:UB);
        tmp{repeat} = reshape(params,dim,num_pts);
        par_matrix(:,:,repeat) = tmp{repeat}';
    end
    %%
    f2 = figure('units','normalized','outerposition',[0 0 1 1]);
    ha = tight_subplot(6,4,[.01 .01],[.1 .01],[.01 .01]);
    param_labels = {'c_3','c_1','a20_{ikk inactivation}','a20_{tnfr inactivation}','k_b','k_f','k_a','k_4','k_{i1}','k_{e1}','k_{t1a}','k_{t2a}','k_{a1a}','k_{d1a}','k_{i2}','k_{e2}','k_{e2a}','c_{4a}','c_{5a}','c_{1a}'};
    %plot figure    
    protein_labels = {'NFkB','IKK','Receptor'};    
    %plot protein parameter distributions
    for i = 1:3
        axes(ha(i));
        p = par_matrix(:,i,:);
        h = histogram(p(:),100,'EdgeColor','none','FaceAlpha',0.5,'FaceColor',[0.8,0.8,0.8]*clr_fr(tcounter)); hold on
        text(0.1,0.8,protein_labels{i},'horizontalAlignment','left','FontSize',fs,'FontName','Arial','Interpreter','none','Units','Normalized');
        set(gca,'FontSize',fs,'FontName','Arial','YTickLabels',[],'XTickLabels',[]);
        hold on
        grid on
        box off
        ylim([0,max(h.Values)]);
        yy = ylim;
        x1 = h.BinEdges((h.Values == max(h.Values)));
        if i <= length(L_par) &&  ~isnan(L_par(i)) 
            xx = [L_par(i),L_par(i)];
            hold on;
            plot(xx,yy,'--r','LineWidth',1.5)
        end    
        if i <= length(R_par) &&  ~isnan(R_par(i)) 
            xx = [R_par(i),R_par(i)];
            hold on;
            plot(xx,yy,'--k','LineWidth',1.5)
        end
        xlim([data.cfg.param_defs{i}.min,data.cfg.param_defs{i}.max])
    end
    axes(ha(4))
    set(gca,'xcolor',[1,1,1],'ycolor',[1,1,1]);
    c = 4;
    subplot_counter = 1;
    umult = unique(multiplier);
    counter=1;
    clr = {[1,0,0],[0,0,1],[0,1,0]};
    for i = 1:length(umult) % for each group, here length(umult) is 3  
        index = find(multiplier==umult(i)); % find in order all parameter indices belonging to a particular group
        for j = index % for each parameter
            par_index = j+3; % add 3 to the index because we have not accounted for the proteins in the multiplier array
            c = c+1;
            %subplot(6,4,c);
            axes(ha(c))
            hold on
            p1 = par_matrix(:,par_index,:); % get the parameter index
            p2 = par_matrix(:,multiplier(j),:); % get the appropriate multiplier 
            p = p1+p2; % add to get compelte parameter value
            h = histogram(p(:),1000,'EdgeColor','none','FaceColor',clr{i}*clr_fr(tcounter),'FaceAlpha',0.5); hold on;
            h.BinEdges(end)
            hold on
            ylim([0 max(h.Values)]);
            yy = ylim;
            % if present, plot the literature estimates
            if par_index <= length(L_par) &&  ~isnan(L_par(par_index)) 
                yy = ylim;
                xx = [L_par(par_index),L_par(par_index)];
                hold on;
                plot(xx,yy,'--r','LineWidth',1.5)
            end
            if par_index <= length(R_par) &&  ~isnan(R_par(par_index)) 
                xx = [R_par(par_index),R_par(par_index)];
                hold on;
                plot(xx,yy,'--k','LineWidth',1.5)
            end
            if par_index <= length(W_par) &&  ~isnan(W_par(par_index)) 
                xx = [W_par(par_index),W_par(par_index)];
                hold on;
                plot(xx,yy,'--b','LineWidth',1.5)
            end           
            xx = [par_min(par_index)+par_min(multiplier(j)),par_min(par_index)+par_min(multiplier(j))];
            hold on;
            plot(xx,yy,'m','LineWidth',1.5)
            xlim([-40 16]); % the x axis range = (LB+LB1, UB+UB1)  = (-35-5,10+6)

            text(0.1,0.8,param_labels{counter},'FontWeight','normal','horizontalAlignment','left','FontSize',fs,'FontName','Arial','Interpreter','Tex','Units','Normalized');
            counter=counter+1;
            if c<=20
                set(gca,'FontSize',fs,'FontName','Arial','YTickLabels',[],'XTickLabels',[]);
            else
                set(gca,'FontSize',fs,'FontName','Arial','YTickLabels',[]);
            end
        end
    end
    %toPDF(f2,['pdf_figures/updated_axisrange_withpublishedpar_NFkB_supplement_posterior_distribution',num2str(trajectory_number),'.pdf'])
end
%%
%fig_print(f2,['NFkB_supplement_posterior_distribution',num2str(trajectory_number),'.png'],[6000,4000],1,1,'w');
%fig_print(f2,'NFkB_supplement_posterior_distribution.png',[6000,4000],1,1,'w');
