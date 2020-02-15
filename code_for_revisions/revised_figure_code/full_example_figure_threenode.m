example_index=1;
addpath('/shared2/LabUserFiles/Sanjana_Gupta/LassoManuscript/lib')
save_pth = '../RevisedFiguresMATLAB/';
pth = '/shared2/LabUserFiles/Sanjana_Gupta/LassoManuscript';
true_params = [-1,-10,0,-10,-10,-10 ] ;
param_labels = {'K_{AB}','K_{BA}','K_{BC}','K_{CB}','K_{CA}','K_{AC}','K_{AB}','K_{BA}','K_{BC}','K_{CB}','K_{CA}','K_{AC}'};
fs = 18;
%%
t = linspace(0,50,8)';
%%
example = 'three_node';
addpath([pth,'/',example])
mu = -10;
B =1;
%load PTLasso data
load([pth,'/',example,'/',example,'_withlasso_mu_',num2str(mu),'_b_',num2str(B),'_repeat_1_progress900000.mat']);
%%
cfg.data{1}.time = t;
cfg.adapt_last = 5e5;
LB = 5e5;
param_indices = 1:length(true_params);
frills = 1;
xlimits = [];
nbins = 1000;
lim = 1;
true_param = true_params;
adjust_y_axis = 0;% 1;
rl_index = find(true_param~=-10);
%Plot parameter posterior distributions
[f1,test_pt_no] = plot_post_separated(cfg,params_chain,energy_chain,param_indices,nbins,true_param,lim,LB,xlimits,frills,rl_index,adjust_y_axis);
%% Plot fits, number of fits = length(LB:step:UB), points per fit = 100
sim_handle1 = str2func(['gensim_',example]);   
n = 100;
step = 100;
sim_handle2 = @(dummy1,dummy2,params)sim_handle1(n,params);
fsim1 = plot_sim(cfg,params_chain,energy_chain,sim_handle2,step,true_param,linspace(0,50,n)',0);
set(gca,'XTickLabels',[],'YTickLabels',[]);
ylim([0 20])
savefig([save_pth,'Figure_threenode_lasso_fit'],'png','-c0','-crop')

set(gca,'FontSize',fs,'FontName','Arial');

%Repeat for without lasso
load([pth,'/',example,'/',example,'_wo_lasso_repeat_1_progress900000.mat']);
%%
cfg.data{1}.time = t;
cfg.adapt_last = 5e5;
LB = 5e5;
f2 = plot_post_separated(cfg,params_chain,energy_chain,param_indices,nbins,true_param,lim,LB,xlimits,frills,rl_index,adjust_y_axis);
fsim2 = plot_sim(cfg,params_chain,energy_chain,sim_handle2,step,true_param,linspace(0,50,n)',0);
ylim([0 20])
set(gca,'XTickLabels',[],'YTickLabels',[]);
savefig([save_pth,'Figure_threenode_wo_lasso_fit'],'png','-c0','-crop')
set(gca,'FontSize',fs,'FontName','Arial');
%%
f_array = [f2,f1];
nplots = 6;

% merge figures into one figure
[hfinal] = merge_posterior_figures_fromhandles(f_array,nplots,[2,6]); 
[hfinal2] = merge_posterior_figures_fromhandles([fsim2,fsim1],1,[2,1]);

h3.fig = figure('units','normalized','outerposition',[0 0 1 1]);
h3.ax = gobjects(15);
h3.ax(2) = subplot(3,6,3:4);
h3.ax(3) = subplot(3,6,5:6);
c = 7;
for iii = 4:15
    h3.ax(iii) = subplot(3,6,c);
    c=c+1;
end
figure_ax = findall(fsim2,'type','axes');   
h3.ax2(2) = copyobj(figure_ax, h3.fig);
h3.ax2(2).Position = h3.ax(2).Position;
ylim([0 20])
delete(h3.ax(2));

figure_ax = findall(fsim1,'type','axes');   
h3.ax2(3) = copyobj(figure_ax, h3.fig);
h3.ax2(3).Position = h3.ax(3).Position;
ylim([0 20])
delete(h3.ax(3));
c = 4;
for iii = 1:12
    h3.ax2(c) = copyobj(hfinal.ax2(iii,1), h3.fig);
    h3.ax2(c).Position = h3.ax(c).Position;
    delete(h3.ax(c));
    c=c+1;
end       
c=1;

h3.ax(2) = subplot(3,6,3:4);
text(0.5,0.8,'Without Lasso','FontSize',fs,'FontName','Arial','Units','normalized','HorizontalAlignment','center');
xlabel('Time (s)');
ylabel({'No. of';'B Molecules'})
set(gca,'FontSize',fs,'FontName','Arial');
h3.ax(3) = subplot(3,6,5:6);
text(0.5,0.8,'With Lasso','FontSize',fs,'FontName','Arial','Units','normalized','HorizontalAlignment','center');
xlabel('Time (s)');
set(gca,'FontSize',fs,'FontName','Arial');
for iii = 7:18
    h3.ax(iii) = subplot(3,6,iii);
    if iii == 7 || iii == 13
        ylabel('Frequency')
    end  
    YY = ylim;
    XX = xlim;
    text(0.5,0.9,param_labels{c},'Color','r','FontSize',fs,'FontName','Arial','Units','normalized','HorizontalAlignment','center');
    set(gca,'FontSize',fs,'FontName','Arial');
    %title(param_labels{c},'FontSize',fs,'FontName','Arial','FontWeight','normal');
    c=c+1;
end   
tightfig(h3.fig);
%savefig('Figure_threenode','png','-c0','-crop')
%toPDF(h3.fig,[save_pth,'/Figure_threenode_correctedJuly23_2019.pdf'])
% fig_print(f1, figname1,[1000,1000],1,1,'w');
