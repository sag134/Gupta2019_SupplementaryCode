addpath('../lib')
addpath('../lib/aboxplot');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
%%
save_pth = './pdf_updated_figures_July23_2019/';
B = [0.1,0.2,0.5,1,5,0.01,0.05];
n = 7e5;
fontsize = 24;
step = 100;
% no of samples used is (N-n:step:N) or n/step, N = length of energy chain
example = 'five_node';
addpath(['../',example])
file_prefix = ['../',example,'/',example,'_withlasso_'];
file_suffix = '_repeat_1_progress900000.mat';
mu_range = [-10,-8];
b_range = sort(B); 
tic   
[mean_gf,std_gf] = get_gf_list(mu_range,b_range,file_prefix,file_suffix,n,step);
toc
f_hyper = hyperparameter_tuning_plot(mean_gf,std_gf,mu_range,b_range);
%savefig('FigureS2_hyperparameter_tuning','png','-c0','-crop');
%toPDF(f_hyper,[save_pth,'/FigureS2_hyperparameter_tuning.pdf'])

% Notes: mu = -10 has intermediate fits that are almost 0 but slowly rise
% and get the last time point; vs. very bad fits that are just flat out 0
% mu = -8, b = 0.5 shows bimodal distributions with good fit and
% intermediate fit. mu = 0.8, b = 0.1 shows bimodal distributions with
% intermediate fit and bad fit
%% get covariation plot for PTLasso
example = 'five_node';
addpath(['../',example])
mu = -10;
true_params =[-1,-10,-10,-10,-10,-10,0,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10];
load(['../',example,'/',example,'_withlasso_mu_',num2str(mu),'_b_',num2str(1),'_repeat_2_progress900000.mat']);
cfg.adapt_last = 2e5;
LB = 2e5;
fcorr = figure;
hold on
scatter3(-10,-10,0,100,'k','filled');
scatter3(0,-10,-10,100,'k','filled');
scatter3(-10,0,-10,100,'k','filled');
p1 = params_chain(1,7,LB:end);
p2 = params_chain(1,9,LB:end);
p3 = params_chain(1,15,LB:end);
xlim([min(p1),max(p1)])
ylim([min(p2),max(p2)])
zlim([min(p3),max(p3)])

xlabel('K_{BC}');
ylabel('K_{BD}');
zlabel('K_{BE}');
set(gca,'FontSize',fontsize,'FontName','Arial');
view(gca,-37.5+180,30)
%toPDF(fcorr,[save_pth,'/FigureS2_correlation_lasso_layer1.pdf'])
hold on
p1 = params_chain(1,7,LB:end);
p2 = params_chain(1,9,LB:end);
p3 = params_chain(1,15,LB:end);
xlim([min(p1),max(p1)])
ylim([min(p2),max(p2)])
zlim([min(p3),max(p3)])
s = scatter3(p1,p2,p3,5,'filled');
s.MarkerFaceColor = [0 0.4470 0.7410];
s.MarkerFaceAlpha = 0.2;  
hold on
scatter3(-10,-10,0,100,'k','filled');
scatter3(0,-10,-10,100,'k','filled');
scatter3(-10,0,-10,100,'k','filled');
set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[],'xlabel',[],'ylabel',[],'zlabel',[])
savefig([save_pth,'/FigureS2_correlation_lasso'],'png','-c0','-crop');
%% Get raw data
load(['../',example,'/',example,'_rawdata.mat']);
xlabel('Time (s)');
ylabel('No. of B molecules');
%savefig('FigureS2_rawdata','png','-c0','-crop');
%toPDF(f,[save_pth,'/FigureS2_rawdata.pdf'])
%% Get covariation plot without lasso
fontsize=24;
load(['../',example,'/',example,'_wo_lasso_repeat_1_progress900000.mat']);
cfg.adapt_last = 2e5;
LB = 2e5;
fcorr = figure;
hold on
scatter3(-10,-10,0,100,'k','filled');
scatter3(0,-10,-10,100,'k','filled');
scatter3(-10,0,-10,100,'k','filled');
p1 = params_chain(1,7,LB:end);
p2 = params_chain(1,9,LB:end);
p3 = params_chain(1,15,LB:end);
xlim([min(p1),max(p1)])
ylim([min(p2),max(p2)])
zlim([min(p3),max(p3)])

xlabel('K_{BC}');
ylabel('K_{BD}');
zlabel('K_{BE}');
set(gca,'FontSize',fontsize,'FontName','Arial');
view(gca,-37.5+180,30)
%toPDF(fcorr,[save_pth,'/FigureS2_correlation_wolasso_layer1.pdf'])
hold on
p1 = params_chain(1,7,LB:end);
p2 = params_chain(1,9,LB:end);
p3 = params_chain(1,15,LB:end);
xlim([min(p1),max(p1)])
ylim([min(p2),max(p2)])
zlim([min(p3),max(p3)])
s = scatter3(p1,p2,p3,5,'filled');
s.MarkerFaceColor = [0 0.4470 0.7410];
s.MarkerFaceAlpha = 0.2;  
hold on
scatter3(-10,-10,0,100,'k','filled');
scatter3(0,-10,-10,100,'k','filled');
scatter3(-10,0,-10,100,'k','filled');
set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[],'xlabel',[],'ylabel',[],'zlabel',[])
savefig([save_pth,'/FigureS2_correlation_wolasso'],'png','-c0','-crop');


%% error comparison box plots
addpath('../five_node');
LB = 2e5;
step = 100;
UB = 9e5;
load('five_node_wo_lasso_repeat_1_progress900000.mat');
e1 = energy_chain(1,LB:step:UB);
load('five_node_withlasso_mu_-10_b_1_repeat_1_progress900000.mat');
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
%savefig('FigureS2_likelihood_comparison','png','-c0','-crop');
%toPDF(f,[save_pth,'/FigureS2_likelihood_comparison.pdf'])