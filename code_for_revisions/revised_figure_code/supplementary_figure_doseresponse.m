param_labels = {'k_{S-RS}','k_{XR-X}','k_{S-XS}','k_{X-0}','k_{R-0}','k_{0-R}'};
addpath('../lib')
addpath('../lib/aboxplot');

% path to parallel tempering scripts
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
% path to supplementary distributions
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/LassoManuscript/lib')
save_pth = '../RevisedFiguresMATLAB/';
pth = '/shared2/LabUserFiles/Sanjana_Gupta/LassoManuscript';

addpath([pth,'/perfect_adaptation_fixedstart/']);
addpath([pth,'/perfect_adaptation_fixedstart/perfect_adaptation/']);
%% PA fits
load([pth,'/perfect_adaptation_fixedstart/perfect_adaptation_continue_chain_two/perfect_adaptation_fixedstart_withlasso_mu_-10_b_1_repeat_1_progress900000.mat'])
true_param = [1,1,0,0,-10,-10,1];
fontsize = 18;
f = figure;
n = 1000;
TT = linspace(0,5,n)';
for d = 1
    for j = 1e5:100:9e5
        params = params_chain(1,:,j);
        [err, sp, obsv] = gensim_perfect_adaptation(n,params);
        p = semilogy([TT;TT(2:end)+TT(end)]',obsv,'b');
        p.Color(4) = 0.01;
        hold on
    end
    T = linspace(0,5,10)';
    errorbar([T;T(2:end)+T(end)]',cfg.data{d}.mean,cfg.data{d}.stdev,'.k','LineWidth',2);
    hold on
    scatter([T;T(2:end)+T(end)]',cfg.data{d}.mean,20,'filled','k')
    [err, sp, obsv] = gensim_perfect_adaptation(n,true_param);
    
    semilogy([TT;TT(2:end)+TT(end)]',obsv,'r','LineWidth',2);
    ylim([0.5 1000])
end
%ylabel('Response (R)');
%xlabel('Time');
%text(0.5,0.9,'With Lasso','Color','r','FontSize',24,'FontName','Arial','Units','normalized','HorizontalAlignment','center');
set(gca,'FontSize',20,'FontName','Arial','XTickLabel',[],'YTickLabel',[]);

savefig([save_pth,'Figure_pa_dr_lasso_fit_axisadjustedDec14'],'png','-c0','-crop')
%fig_print(f,'perfect_adaptation_fit_lasso_S3B.png',[1000,1000],1,1,'w');
%%toPDF(f,'pdf_figures/perfect_adaptation_fit_lasso_S3B.pdf');
%%
load([pth,'/perfect_adaptation_fixedstart/perfect_adaptation_continue_chain_two/perfect_adaptation_fixedstart_wo_lasso_repeat_1_progress900000.mat'])
true_param = [1,1,0,0,-10,-10,1];
f = figure;
n = 1000;
TT = linspace(0,5,n)';
t = linspace(0,800,n)';
for d = 1
    for j = 1e5:100:9e5
        params = params_chain(1,:,j);
        [err, sp, obsv] = gensim_perfect_adaptation(n,params);
        p = semilogy([TT;TT(2:end)+TT(end)]',obsv,'b');
        p.Color(4) = 0.01;
        hold on
    end
    T = linspace(0,5,10)';
    errorbar([T;T(2:end)+T(end)]',cfg.data{d}.mean,cfg.data{d}.stdev,'.k','LineWidth',2);
    hold on
    scatter([T;T(2:end)+T(end)]',cfg.data{d}.mean,20,'filled','k')
    [err, sp, obsv] = gensim_perfect_adaptation(n,true_param);   
    semilogy([TT;TT(2:end)+TT(end)]',obsv,'r','LineWidth',2);
    ylim([0.5 1000])
end
%ylabel('Response (R)');
%xlabel('Time');
%text(0.5,0.9,'Without Lasso','Color','r','FontSize',24,'FontName','Arial','Units','normalized','HorizontalAlignment','center');
set(gca,'FontSize',20,'FontName','Arial');
set(gca,'FontSize',20,'FontName','Arial','XTickLabel',[],'YTickLabel',[]);
savefig([save_pth,'Figure_pa_dr_wolasso_fit_Axisadjusted_Dec14'],'png','-c0','-crop')
%fig_print(f,'perfect_adaptation_fit_wolasso_S3A.png',[1000,1000],1,1,'w');
%toPDF(f,'pdf_figures/perfect_adaptation_fit_wolasso_S3A.pdf');
%% Covariance plot
load([pth,'/perfect_adaptation_fixedstart/perfect_adaptation_continue_chain_two/perfect_adaptation_fixedstart_withlasso_mu_-10_b_1_repeat_1_progress900000.mat'])
p1 = params_chain(1,2,1e5:9e5);
p2 = params_chain(1,3,1e5:9e5);
f = figure;scatter(p1,p2,20,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.005,'MarkerEdgeColor','none');
xlim([-4 6]);
ylim([-6 5]);
set(gca,'FontSize',fontsize);
xlabel(param_labels{2});
ylabel(param_labels{3});
%fig_print(f,'perfect_adaptation_covariance_S3C.png',[1000,1000],1,1,'w');
%%toPDF(f,'pdf_figures/perfect_adaptation_covariance_S3C.pdf');
%% Hyperparameter tuning perfect adaptation
n = 8e5;
file_prefix = [pth,'/perfect_adaptation_fixedstart/perfect_adaptation_continue_chain_two/perfect_adaptation_fixedstart_withlasso_'];
file_suffix = '_repeat_1_progress900000.mat';
mu_range = [-10,-8];
b_range =  sort([5,1,0.5,0.1]);%sort(B{i});
tic
[mean_gf,std_gf] = get_gf_list(mu_range,b_range,file_prefix,file_suffix,n,1000);
toc
f = hyperparameter_tuning_plot(mean_gf,std_gf,mu_range,b_range);
figname = ['perfect_adaptation_hyperparameter_tuning.png'];
%fig_print(f, figname,[1000,1000],1,1,'w');
%toPDF(f,'pdf_updated_figures_July23_2019/perfect_adaptation_hyperparameter_tuning.pdf');
%% Hyperparameter linear dose response

% Note: In this example regularization with b >=0.5 gives good fits,
% regularization with b >0.1 gives bad fits (all parameters are 0 so model output is 0). for b = 0.1 we get an
% intermediate fit where the reaction S->S+R is the only reaction giving a
% linear model output for R that fits the data poorly.

n = 4e5;
file_prefix = [pth,'/linear_full_trajectory/linear_full_trajectory_withlasso_'];
addpath([pth,'/linear_full_trajectory/']);
file_suffix = '_repeat_1_progress900000.mat';
mu_range = [-10,-8];
b_range =  sort([1,0.5,5,0.01,0.05,0.1]);%sort(B{i});
tic
[mean_gf,std_gf] = get_gf_list(mu_range,b_range,file_prefix,file_suffix,n,1000);
toc
f = hyperparameter_tuning_plot(mean_gf,std_gf,mu_range,b_range);
figname = ['linear_hyperparameter_tuning.png'];
%fig_print(f, figname,[1000,1000],1,1,'w');
%toPDF(f,'pdf_updated_figures_July23_2019/linear_hyperparameter_tuning.png.pdf');
%%
addpath([pth,'/linear_full_trajectory/']);
step = 1000;
LB = 5e5;
UB = 9e5;
load('linear_full_trajectory_wo_lasso_repeat_1_progress900000.mat');
e1 = energy_chain(1,LB:step:UB);
load('linear_full_trajectory_withlasso_mu_-10_b_0.5_repeat_1_progress900000.mat');
e2 = zeros(1,length(LB:step:UB));
c = 1;
for i = LB:step:UB
    i
    params = params_chain(1,:,i);
    e2(c) = likelihood(params,cfg);
    c=c+1;
end
f = figure;
hold on
aboxplot({e1',e2'},'OutlierMarker','none')

%Significance testing not used for paper.
[p,h] = ranksum(e1,e2)
if h==1 && p<1e-4
    [q1 q2 q3 fu1 fl1 ou ol] = quartile(e1);
    [q1 q2 q3 fu2 fl2 ou ol] = quartile(e2);
    hold on
    hold on
    plot([0.825,1.175],[fu1,fu1]+2,'LineWidth',2,'Color',[0.5,0.5,0.5])
    hold on
    plot([0.825,0.825],[fu1+0.2,fu1+2],'LineWidth',2,'Color',[0.5,0.5,0.5])
    plot([1.175,1.175],[fu2+0.2,fu1+2],'LineWidth',2,'Color',[0.5,0.5,0.5])
    hold on
    scatter(0.825+(1.175-0.825)/2,fu1+3,150,'*','MarkerEdgeColor',[0.3,0.3,0.3])
end
ylim([0 8])
%toPDF(f,'pdf_updated_figures_July23_2019/FigureS3_linear_likelihood_comparison.pdf')
%%
addpath('../perfect_adaptation_fixedstart/perfect_adaptation_continue_chain_two/');
step = 1000;
LB = 1e5;
UB = 9e5;
load('perfect_adaptation_fixedstart_wo_lasso_repeat_1_progress900000.mat');
e1 = energy_chain(1,LB:step:UB);
load('perfect_adaptation_fixedstart_withlasso_mu_-10_b_1_repeat_1_progress900000.mat');
e2 = zeros(1,length(LB:step:UB));
c = 1;
for i = LB:step:UB
    i
    params = params_chain(1,:,i);
    e2(c) = likelihood(params,cfg);
    c=c+1;
end

f = figure;
hold on
aboxplot({e1',e2'},'OutlierMarker','none')

%Significance testing not used for paper.
[p,h] = ranksum(e1,e2)    
if h==1 && p<1e-4
    [q1 q2 q3 fu1 fl1 ou ol] = quartile(e1);
    [q1 q2 q3 fu2 fl2 ou ol] = quartile(e2);
    hold on
    hold on
    plot([0.825,1.175],[fu1,fu1]+2,'LineWidth',2,'Color',[0.5,0.5,0.5])
    hold on
    plot([0.825,0.825],[fu1+0.2,fu1+2],'LineWidth',2,'Color',[0.5,0.5,0.5])
    plot([1.175,1.175],[fu2+0.2,fu1+2],'LineWidth',2,'Color',[0.5,0.5,0.5])
    hold on
    scatter(0.825+(1.175-0.825)/2,fu1+3,150,'*','MarkerEdgeColor',[0.3,0.3,0.3])
end
ylim([0 10])
%toPDF(f,'pdf_updated_figures_July23_2019/FigureS3_PA_likelihood_comparison.pdf')
