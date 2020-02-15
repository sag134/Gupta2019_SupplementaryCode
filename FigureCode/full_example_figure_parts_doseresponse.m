addpath('../lib')
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
param_labels = {'k_{S-RS}','k_{XR-X}','k_{S-XS}','k_{X-0}','k_{R-0}','k_{0-R}'};
true_param = [1,-10,-10,-10,-2,-10,-10];
fontsize = 18;
%%
addpath('../linear_full_trajectory');
load('../linear_full_trajectory/linear_full_trajectory_withlasso_mu_-10_b_0.5_repeat_1_progress900000.mat')
%%
param_indices = 1:6;
cfg.adapt_last = 5e5;
LB = 5e5;
frills = 1;
xlimits = [];
nbins = 1000;
lim = 1;
rl_index = find(true_param~=-10);
adjust_y_axis = 0;
f1 = plot_post_separated(cfg,params_chain,energy_chain,param_indices,nbins,true_param,lim,LB,xlimits,frills,rl_index,adjust_y_axis);
figure_ax = findall(f1,'type','axes');
h3.fig = figure('units','normalized','outerposition',[0 0 1 1]);
c = 1;
for j = 6:-1:1
    h3.ax(j) = subplot(1,6,j); 
    h3.ax2(j) = copyobj(figure_ax(c,1), h3.fig);
    h3.ax2(j).Position = h3.ax(j).Position;
    delete(h3.ax(j));
    text(0.5,0.8,param_labels{j},'Color','r','FontSize',fontsize,'FontName','Arial','Units','normalized','HorizontalAlignment','center');
    c=c+1;
    if j == 1
        ylabel('Frequency');
    end
    set(gca,'FontSize',22,'FontName','Arial');
end
%%
%fig_print(h3.fig,'linear_posterior_lasso_Figure3D_bottom.png',[5000,1000],1,1,'w');
%toPDF(h3.fig,'pdf_figures/linear_posterior_lasso_Figure3D_bottom.pdf');

f2 = figure;
n = 1000;
t = linspace(0,800,n)';
for d = 1:4
    for j =  LB:1000:cfg.nswaps
        params = params_chain(1,:,j);
        [err, sp, obsv] = gensim_linear(d,n,params);
        p = plot(t,obsv,'b');
        p.Color(4) = 0.1;
        hold on
    end
    errorbar(linspace(0,800,4)',cfg.data{d}.mean,cfg.data{d}.stdev,'.k','LineWidth',2);
    scatter(linspace(0,800,4)',cfg.data{d}.mean,20,'filled','k')
    [err, sp, obsv] = gensim_linear(d,n,true_param);
    plot(t,obsv,'r','LineWidth',2);
end
ylabel('Response (R)');
xlabel('Time');
text(0.5,0.9,'With Lasso','Color','r','FontSize',24,'FontName','Arial','Units','normalized','HorizontalAlignment','center');
set(gca,'FontSize',20,'FontName','Arial');
figname2 = 'linear_fit_lasso_FigureC.png';
%fig_print(f2, figname2,[5000,2000],1,1,'w');
%toPDF(f2,'pdf_figures/linear_fit_lasso_FigureC.pdf');
%%
load('../linear_full_trajectory/linear_full_trajectory_wo_lasso_repeat_1_progress900000.mat')
param_indices = 1:6;
cfg.adapt_last = 5e5;
LB = 5e5;
frills = 1;
xlimits = [];
nbins = 1000;
lim = 1;
true_param = true_params{3};
rl_index = find(true_param~=-10);
adjust_y_axis = 0;
f1 = plot_post_separated(cfg,params_chain,energy_chain,param_indices,nbins,true_param,lim,LB,xlimits,frills,rl_index,adjust_y_axis);
figure_ax = findall(f1,'type','axes');
h3.fig = figure('units','normalized','outerposition',[0 0 1 1]);
c = 1;
for j = 6:-1:1
    h3.ax(j) = subplot(1,6,j); 
    h3.ax2(j) = copyobj(figure_ax(c,1), h3.fig);
    h3.ax2(j).Position = h3.ax(j).Position;
    delete(h3.ax(j));
    text(0.5,0.8,param_labels{j},'Color','r','FontSize',fontsize,'FontName','Arial','Units','normalized','HorizontalAlignment','center');
    c=c+1;
    if j == 1
        ylabel('Frequency');
    end
    set(gca,'FontSize',22,'FontName','Arial');
end
%fig_print(h3.fig,'linear_posterior_wolasso_Figure3D_top.png',[5000,1000],1,1,'w');
%toPDF(h3.fig,'pdf_figures/linear_posterior_wolasso_Figure3D_top.pdf');

f = figure;
n = 1000;
t = linspace(0,800,n)';
for d = 1:4
    for j = LB:1000:cfg.nswaps
        params = params_chain(1,:,j);
        [err, sp, obsv] = gensim_linear(d,n,params);
        p = plot(t,obsv,'b');
        p.Color(4) = 0.1;
        hold on
    end
    errorbar(linspace(0,800,4)',cfg.data{d}.mean,cfg.data{d}.stdev,'.k','LineWidth',2);
    scatter(linspace(0,800,4)',cfg.data{d}.mean,20,'filled','k')
    [err, sp, obsv] = gensim_linear(d,n,true_param);
    plot(t,obsv,'r','LineWidth',2);
end
ylabel('Response (R)');
xlabel('Time');
text(0.5,0.9,'Without Lasso','Color','r','FontSize',24,'FontName','Arial','Units','normalized','HorizontalAlignment','center');
set(gca,'FontSize',20,'FontName','Arial');
figname2 = 'linear_fit_wolasso_FigureB.png';
%fig_print(f, figname2,[5000,2000],1,1,'w');
%toPDF(f,'pdf_figures/linear_fit_wolasso_FigureB.pdf');
%%
addpath('../perfect_adaptation_fixedstart');
load('../perfect_adaptation_fixedstart/perfect_adaptation_continue_chain_two/perfect_adaptation_fixedstart_withlasso_mu_-10_b_1_repeat_1_progress900000.mat')
param_indices = 1:6;
cfg.adapt_last = 1e5;
LB = 1e5;
frills = 1;
xlimits = [];
nbins = 1000;
lim = 1;
true_param = true_params{4};
rl_index = find(true_param~=-10);
adjust_y_axis = 0;
f2 = plot_post_separated(cfg,params_chain,energy_chain,param_indices,nbins,true_param,lim,LB,xlimits,frills,rl_index,adjust_y_axis);
%%
figure_ax = findall(f2,'type','axes');
h3.fig = figure('units','normalized','outerposition',[0 0 1 1]);
c = 1;
for j = 6:-1:1
    h3.ax(j) = subplot(1,6,j); 
    h3.ax2(j) = copyobj(figure_ax(c,1), h3.fig);
    h3.ax2(j).Position = h3.ax(j).Position;
    delete(h3.ax(j));
    text(0.5,0.8,param_labels{j},'Color','r','FontSize',fontsize,'FontName','Arial','Units','normalized','HorizontalAlignment','center');
    c=c+1;
    if j == 1
        ylabel('Frequency');
    end
    set(gca,'FontSize',22,'FontName','Arial');
end

%fig_print(h3.fig,'perfect_adaptation_posterior_lasso_Figure3F_bottom.png',[5000,1000],1,1,'w');
%toPDF(h3.fig,'pdf_figures/perfect_adaptation_posterior_lasso_Figure3F_bottom.pdf');

%%
load('../perfect_adaptation_fixedstart/perfect_adaptation_continue_chain_two/perfect_adaptation_fixedstart_wo_lasso_repeat_1_progress900000.mat')
param_indices = 1:6;
cfg.adapt_last = 1e5;
LB = 1e5;
frills = 1;
xlimits = [];
nbins = 1000;
lim = 1;
true_param = true_params{4};
rl_index = find(true_param~=-10);
adjust_y_axis = 0;
f2 = plot_post_separated(cfg,params_chain,energy_chain,param_indices,nbins,true_param,lim,LB,xlimits,frills,rl_index,adjust_y_axis);

figure_ax = findall(f2,'type','axes');
h3.fig = figure('units','normalized','outerposition',[0 0 1 1]);
c = 1;
for j = 6:-1:1
    h3.ax(j) = subplot(1,6,j); 
    h3.ax2(j) = copyobj(figure_ax(c,1), h3.fig);
    h3.ax2(j).Position = h3.ax(j).Position;
    delete(h3.ax(j));
    text(0.5,0.8,param_labels{j},'Color','r','FontSize',fontsize,'FontName','Arial','Units','normalized','HorizontalAlignment','center');
    c=c+1;
    if j == 1
        ylabel('Frequency');
    end
    set(gca,'FontSize',22,'FontName','Arial');
end
%fig_print(h3.fig,'perfect_adaptation_posterior_WOlasso_Figure3F_top.png',[5000,1000],1,1,'w');
toPDF(h3.fig,'pdf_figures/perfect_adaptation_posterior_WOlasso_Figure3F_top.pdf');

% Note: PA fit figure code is in the supplementary figure code