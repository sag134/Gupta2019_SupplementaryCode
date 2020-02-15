example_index=1;
addpath('/shared2/LabUserFiles/Sanjana_Gupta/LassoManuscript/lib')
save_pth = '../RevisedFiguresMATLAB/';
true_params = [-1,-10,0,-10,-10,-10 ] ;
%%
param_labels = {'K_{AB}','K_{BA}','K_{BC}','K_{CB}','K_{CA}','K_{AC}','K_{AB}','K_{BA}','K_{BC}','K_{CB}','K_{CA}','K_{AC}'};
fs = 18;
pth = '../three_node_parameter_variation';
addpath(pth);
data = load('../three_node_parameter_variation/full_three_node_parameter_variability0.05.mat');
parameters = [ -1, -10, 0, -10, -10, -10 ]; 
f = figure('Units','normalized','OuterPosition',[0,0,1,0.5]);
c=1;
subplot(1,2,c);c=c+1;
p1 = data.params(:,1);
p2 = data.params(:,3);
scatter(p1,p2,50,'filled','k');
hold on
if sum(data.parameters~=parameters)~=0
    disp('PROBLEM')
end
scatter(data.parameters(1),data.parameters(3),'filled','r');


subplot(1,2,c);
t = linspace(0,50,8)';
plot(t,data.data,'LineWidth',1.5);
hold on
[err, sp, obsv] = simulate_three_node([],[],parameters);
hold on
errorbar(t,data.expt{1}.mean,data.expt{1}.stdev,'.k','LineWidth',1.5);
hold on
scatter(t,data.expt{1}.mean,50,'k','filled')
plot(t,obsv,'--r','LineWidth',1.5);
hold on
ylim([0 6.5])
%%toPDF(f,[save_pth,'variable_parameters_data_Dec16.pdf'])

%%
t = linspace(0,50,8)';
pth = '../three_node_parameter_variation/';
addpath(pth);
load([pth,'three_node_parameter_variability_withlasso_fraction_error_5_mu_-10_b_1_repeat_1_progress900000.mat'])
cfg.data{1}.time = t;
%%
cfg.adapt_last = 1e5;
LB = 1e5;
param_indices = 1:length(true_params);
frills = 1;
xlimits = [];
nbins = 1000;
lim = 1;
true_param = true_params;
adjust_y_axis = 0;
rl_index = find(true_param~=-10);
%%
%Plot parameter posterior distributions
[f1,test_pt_no] = plot_post_separated(cfg,params_chain,energy_chain,param_indices,nbins,true_param,lim,LB,xlimits,frills,rl_index,adjust_y_axis);
nplots = 6;
[h] = merge_posterior_figures_fromhandles(f1,nplots,[1,6]); 
%%toPDF(h.fig,[save_pth,'variable_parameter_data_threenodedistributionDec16.pdf'])
%set(gca,'XTickLabels',[],'YTickLabels',[]);
%savefig([save_pth,'variable_parameter_data_threenodedistributionDec16'],'png','-c0','-crop')
%%
sim_handle1 = str2func(['gensim_three_node']);   
n = 100;
step = 100;
sim_handle2 = @(dummy1,dummy2,params)sim_handle1(n,params);
fsim1 = plot_sim(cfg,params_chain,energy_chain,sim_handle2,step,true_param,linspace(0,50,n)',0);
ylim([0 12])
%%
set(gca,'FontSize',fs,'FontName','Arial');
%%toPDF(fsim1,[save_pth,'variable_parameter_data_threenodefitDec16.pdf'])
set(gca,'XTickLabels',[],'YTickLabels',[]);
%savefig([save_pth,'variable_parameter_data_threenodefitDec16'],'png','-c0','-crop')