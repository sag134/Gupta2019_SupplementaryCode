addpath('../SingleCellNFkB_pulse');
%path to parallel tempering files
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
addpath('../lib');
counter=1;
LB = 5e4; %for biorxiv version this was set to 1e4. Changing it to 5e4 for consistency with other figures
UB = 9.9e5;
step = 5e2;%number of samples per repeat. Total number of parameters sampled here = 5e2*
mu = -25;
b = 2;
repeat = 1;
trajectory_number = 2; %a representative trajectory
tequil = (0:1e5:1e8)';
tstim = (0:1:5*60)';
dose = 5;
twash = (0:5*60:(48-1)*5*60)';
t = tstim(end)+twash;
data = load(['../SingleCellNFkB_pulse/SingleCellNFkB_continue_chain_reduced_model_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat),'_progress990000.mat']);    
cfg = data.cfg;
%%
neqpts = 1; %number of equilibration points to plot
figure
f2 = figure('units','normalized','outerposition',[0 0 1 1]);
fs = 12;
for j = round(linspace(LB,UB,step))
    display(['parameter set: ',num2str(j)]);
    parameters =  data.params_chain(1,:,j);
    [err, sp, obsv,oeq,o1,o2] = predict( [],[], parameters,cfg);
    observable_labels = { 'nNFKB', 'cNFKB', 'cikba', 'Ccomplex', 'nikba', 'ncomplex', 'ikkN', 'ikka', 'TNFR_a', 'TNFR_i', 'A20', 'TNF' };
    for i = 1:12
        subplot(3,4,i);
        if i ==11 %Plot A20 trajectories in log scale
            p = semilogy([tstim(1:end-1)/3600;t/(3600)],[o1(1:end-1,i)',o2(:,i)'],'LineWidth',2);
        else
            p = plot([tstim(1:end-1)/3600;t/(3600)],[o1(1:end-1,i)',o2(:,i)'],'LineWidth',2);
        end
        p.Color(4) = 0.5;
        hold on
        set(gca,'XTickLabels',[]);
    end
    hold on
end

for i = 1:12
    subplot(3,4,i);
    yy = ylim;
    yy(1) = 0;
    yy
    ylim(yy);
end
toPDF(f2,'nfkb_updated_figures_Sept3/nfkb_model_prediction_sept2.pdf');