example_index=1;
addpath('/shared2/LabUserFiles/Sanjana_Gupta/LassoManuscript/lib')
save_pth = '../RevisedFiguresMATLAB/';
true_params = [ -1, -20, -20, -20, -20, -20, 0, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20 ];

t = linspace(0,50,8)';
%%
pth = '../five_node_parameter_variation/';
addpath(pth);
load([pth,'five_node_withlasso_mu_-10_b_1_repeat_1_progress900000.mat'])
%%
cfg.data{1}.time = t;
cfg.adapt_last = 1e5;
LB = 1e5;
param_indices = 1:length(true_params);
frills = 1;
xlimits = [];
nbins = 1000;
lim = 1;
true_param = true_params;
adjust_y_axis = 0;
rl_index = find(true_param~=-20);
%Plot parameter posterior distributions
[f1,test_pt_no] = plot_post_separated(cfg,params_chain,energy_chain,param_indices,nbins,true_param,lim,LB,xlimits,frills,rl_index,adjust_y_axis);
nplots = 20;
[h] = merge_posterior_figures_fromhandles(f1,nplots,[2,10]); 
%toPDF(h.fig,[save_pth,'variable_parameter_data_fivenodedistributionDec16.pdf'])
%%
sim_handle1 = str2func(['gensim_five_node']);   
n = 100;
step = 100;
sim_handle2 = @(dummy1,dummy2,params)sim_handle1(n,params);
fsim1 = plot_sim(cfg,params_chain,energy_chain,sim_handle2,step,true_param,linspace(0,50,n)',0);
ylim([0 12])
%set(gca,'FontSize',fs,'FontName','Arial');
%toPDF(fsim1,[save_pth,'variable_parameter_data_fivenodefitDec16.pdf'])
set(gca,'XTickLabels',[],'YTickLabels',[]);
%savefig([save_pth,'variable_parameter_data_fivenodefitDec16'],'png','-c0','-crop')