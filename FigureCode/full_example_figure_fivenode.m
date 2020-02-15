example_index=2;
addpath('../lib')
save_pth = 'pdf_updated_figures_July23_2019';
true_params = [-1,-10,-10,-10,-10,-10,0,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10];
%%
B = 1;
mu = -10;
param_labels = {'K_{AB}','K_{BA}','K_{AC}','K_{CA}','K_{AD}','K_{DA}','K_{BC}','K_{CB}','K_{BD}','K_{DB}','K_{CD}','K_{DC}','K_{AE}','K_{EA}','K_{BE}','K_{EB}','K_{CE}','K_{EC}','K_{DE}','K_{ED}','K_{AB}','K_{BA}','K_{AC}','K_{CA}','K_{AD}','K_{DA}','K_{BC}','K_{CB}','K_{BD}','K_{DB}','K_{CD}','K_{DC}','K_{AE}','K_{EA}','K_{BE}','K_{EB}','K_{CE}','K_{EC}','K_{DE}','K_{ED}'};
fontsize = 18;

example = 'five_node';
addpath(['../',example])
load(['../',example,'/',example,'_withlasso_mu_',num2str(mu),'_b_',num2str(B),'_repeat_1_progress900000.mat']);
%%
cfg.adapt_last = 2e5;
cfg.data{1}.time = linspace(0,50,8);
LB = 2e5;
param_indices = 1:length(true_params);
frills = 1;
xlimits = [];
nbins = 1000;
lim = 1;
true_param = true_params;
adjust_y_axis = 0;
rl_index = find(true_param~=-10);
%plot posterior parameter distribution with PTLasso
[f1,test_pt_no] = plot_post_separated(cfg,params_chain,energy_chain,param_indices,nbins,true_param,lim,LB,xlimits,frills,rl_index,adjust_y_axis);
%% Plot fits
sim_handle1 = str2func(['gensim_',example]);   
n = 100;
step = 100;
sim_handle2 = @(dummy1,dummy2,params)sim_handle1(n,params);
fsim1 = plot_sim(cfg,params_chain,energy_chain,sim_handle2,step,true_param,linspace(0,50,n)',0);
ylim([0 15]);
%% without lasso posterior parameter distribution and fits
load(['../',example,'/',example,'_wo_lasso_repeat_1_progress900000.mat']);
cfg.adapt_last = 2e5;
cfg.data{1}.time = linspace(0,50,8);
LB = 2e5;
f2 = plot_post_separated(cfg,params_chain,energy_chain,param_indices,nbins,true_param,lim,LB,xlimits,frills,rl_index,adjust_y_axis);
fsim2 = plot_sim(cfg,params_chain,energy_chain,sim_handle2,step,true_param,linspace(0,50,n)',0);
ylim([0 80]);
f_array = [f2,f1];
nplots = 20;
[hfinal] = merge_posterior_figures_fromhandles(f_array,nplots,[4,10]);
%%
[hfinal2] = merge_posterior_figures_fromhandles([fsim2,fsim1],1,[2,1]);
%% put figures together
h3.fig = figure('units','normalized','outerposition',[0 0 1 1]);
h3.ax = gobjects(43);
%h3.ax(1) = subplot(3,6,1:2);
h3.ax(2) = subplot(6,10,[4,5,6,14,15,16]);
figure_ax = findall(fsim2,'type','axes');   
h3.ax2(2) = copyobj(figure_ax, h3.fig);
h3.ax2(2).Position = h3.ax(2).Position;
delete(h3.ax(2));
text(0.5,0.8,'Without Lasso','FontSize',fontsize,'FontName','Arial','Units','normalized','HorizontalAlignment','center');

xlabel('Time (s)');
ylabel({'No. of B';'Molecules'})

h3.ax(3) = subplot(6,10,[8,9,10,18,19,20]);
figure_ax = findall(fsim1,'type','axes');   
h3.ax2(3) = copyobj(figure_ax, h3.fig);
h3.ax2(3).Position = h3.ax(3).Position;


delete(h3.ax(3));
text(0.5,0.8,'With Lasso','FontSize',fontsize,'FontName','Arial','Units','normalized','HorizontalAlignment','center');
xlabel('Time (s)');
ylabel({'No. of B';'Molecules'})

c = 1;
c1 = 21;
for iii = 4:43
    h3.ax(iii) = subplot(6,10,c1);
    c1=c1+1;


    h3.ax2(iii) = copyobj(hfinal.ax2(c,1), h3.fig);
    h3.ax2(iii).Position = h3.ax(iii).Position;
    delete(h3.ax(iii));
    h3.ax2(iii);
    text(0.5,0.8,param_labels{c},'Color','r','FontSize',fontsize,'FontName','Arial','Units','normalized','HorizontalAlignment','center');
    if mod((c-1),10)==0
        ylabel('Frequency','FontSize',fontsize','FontName','Arial');
    end

    YY = ylim;
    XX = xlim;
    c=c+1;  
end
%%

tightfig(h3.fig);
toPDF(h3.fig,[save_pth ,'/Figure_fivenode_July23_2019'])
%toPDF(h3.fig,'pdf_figures/Figure_fivenode_fixedsimfigures')
%savefig('Figure_fivenode','png','-c0','-crop')
%fig_print(f1, figname1,[1000,1000],1,1,'w');
