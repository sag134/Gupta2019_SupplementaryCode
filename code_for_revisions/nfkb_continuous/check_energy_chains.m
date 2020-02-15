%Check the energy chains visually to make sure they look flat
f = figure('units','normalized','outerposition',[0,0,1,1]);
c=1;
for i = 1:16
    i
    load(['SingleCellNFkB_reduced_model_continuous_adaptiveIKK_fixedstart_withlasso_mu_-25_b_2_trajectory_1_repeat_',num2str(i),'_progress500000.mat'])
    subplot(8,2,c);c=c+1;
    plot(energy_chain(2,1e5:5e5));hold on
    e = energy_chain(1,1e5:5e5);
    %histogram(e(:),'EdgeColor','none');
    %hold on
    plot(e);
    hold on   
    ylim([20 50])
%return
end
