addpath('../lib')
addpath('../lib/aboxplot/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
save_pth = './';
B = [0.1,0.5,1,5,0.01,0.05];
n = 4e5;
step = 100;
fontsize = 22;
t = linspace(0,50,8)';
folder_name = 'pdf_updated_figures_July23_2019';
example = 'three_node';
addpath(['../',example]) %path to the data files / directory from which PT or PTLasso was run 
%%
file_prefix = ['../',example,'/',example,'_withlasso_'];
file_suffix = '_repeat_1_progress900000.mat';
%%
mu_range = [-10,-8];
b_range = sort(B); 
%hyperparameter tuning figure for n/step PTLasso samples
tic   
[mean_gf,std_gf] = get_gf_list(mu_range,b_range,file_prefix,file_suffix,n,step);
toc
f_hyper = hyperparameter_tuning_plot(mean_gf,std_gf,mu_range,b_range);
%toPDF(f_hyper,[folder_name,'/FigureS1_hyperparameter_tuning.pdf'])

% Note: standard deviation during cross over point seems higher because of
% bimodal distribution, high likelihood low prior and high prior low
% likelihood give the same energy

mu = -10;
true_param = [-1,-10,0,-10,-10,-10];
load(['../',example,'/',example,'_withlasso_mu_',num2str(mu),'_b_',num2str(0.1),'_repeat_1_progress900000.mat']);
cfg.adapt_last = 5e5;
LB = 5e5;
param_indices = 1:length(true_param);
frills = 1;
xlimits = [];
nbins = 1000;
lim = 1;
adjust_y_axis = 0;% 1;
rl_index = find(true_param~=-10);
%plot posteriors
%adjust y axis = 0 means that the code won't go through and make all the y
%axes the same value; Good thing to emphasize shapes instead of heights
[f1,test_pt_no] = plot_post_separated_subplotordering(cfg,params_chain,energy_chain,param_indices,nbins,true_param,lim,LB,xlimits,frills,rl_index,adjust_y_axis,[1 6]);
%toPDF(f1,[folder_name,'S1D_posteriors_lowb.pdf'])
%%

% plot fits, gensim_threenode takes as input (n,params) where n is the step
% number of parameters to sample = length(LB:n:UB)
sim_handle1 = str2func(['gensim_',example]);   
n = 100;
step = 100;
sim_handle2 = @(dummy1,dummy2,params)sim_handle1(n,params);
cfg.data{1}.time = t;
fsim1 = plot_sim(cfg,params_chain,energy_chain,sim_handle2,step,true_param,linspace(0,50,n)',0);
xlabel('Time (s)');
ylabel({'No. of';'B Molecules'})
ylim([0 8])
set(gca,'FontSize',fontsize,'FontName','Arial');
%toPDF(fsim1,[folder_name,'/S1C_lowb_fits.pdf'])
%% plot raw data saved already
load(['../',example,'/',example,'_rawdata.mat']);
%fig_print(f,'S1A_data.png',[1000,1000],1,1,'w');
xlabel('Time (s)');
ylabel('No. of B molecules');
%savefig('FigureS2_rawdata','png','-c0','-crop');
%toPDF(f,[folder_name,'/FigureS1A_rawdata.pdf'])
%% Get comparison of log likelihood distribution
addpath('../three_node');
LB = 5e5;
step = 100;
UB = 9e5;
load('three_node_wo_lasso_repeat_1_progress900000.mat');
e1 = energy_chain(1,LB:step:UB);
load('three_node_withlasso_mu_-10_b_1_repeat_1_progress900000.mat');
e2 = zeros(1,length(LB:step:UB));
c = 1;
for i = LB:step:UB
    i
    params = params_chain(1,:,i);
    e2(c) = likelihood(params,cfg);
    c=c+1;
end
%aboxplot({e1',e2'})
fs  = 8;
%%
f = figure;
set(gca,'XTickLabels',[],'FontSize',8);
hold on
aboxplot({e1',e2'},'outliermarker','none')
ylim([0 5])
[h,p] = ttest2(e1',e2');
[p,h] = ranksum(e1',e2');
% stuff to plot the significance if desired (we don't use this for the paper)
% all the numbers are heuristics I adjusted for plotting purposes (location
% of the significance star etc
if h==1
    [q1 q2 q3 fu1 fl1 ou ol] = quartile(e1);
    [q1 q2 q3 fu2 fl2 ou ol] = quartile(e2);
    hold on
    plot([0.825,1.175],[fu1,fu1]+0.8,'LineWidth',2,'Color',[0.5,0.5,0.5])
    hold on
    plot([0.825,0.825],[fu1+0.2,fu1+0.8],'LineWidth',2,'Color',[0.5,0.5,0.5])
    plot([1.175,1.175],[fu2+0.2,fu1+0.8],'LineWidth',2,'Color',[0.5,0.5,0.5])
    hold on
    scatter(0.825+(1.175-0.825)/2,fu1+1,150,'*','MarkerEdgeColor',[0.3,0.3,0.3])   
end
%toPDF(f,[folder_name,'/FigureS1_likelihood_comparison.pdf'])