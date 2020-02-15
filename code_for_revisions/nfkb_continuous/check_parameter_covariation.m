%% This is not included in the paper but here we are checking to make sure correct parameter covariations were maintained for the continuous treatment runs.
for i = 1:16
    load(['SingleCellNFkB_reduced_model_continuous_adaptiveIKK_fixedstart_withlasso_mu_-25_b_2_trajectory_1_repeat_',num2str(i),'_progress500000.mat'])
    %%
    par_matrix = params_chain(1,:,1e5:5e5);
    f2 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1);

    p1 = par_matrix(:,13,:)+par_matrix(:,26,:); % IkB module
    p2 = par_matrix(:,14,:)+par_matrix(:,26,:); % IkB module

    s = binscatter(p1(:),p2(:),100);
    hold on
    scatter(p1(:),p2(:),1,'filled')
    xlim([-20 10]);ylim([-20 10])

    hold on
    plot([-20 10],[-20,10],'--k','LineWidth',2)
    xlabel('Exit rate of free nuclear IkB')
    ylabel('Exit rate of nuclear IkB-NFkB complex');


    subplot(2,2,2);

    p1 = par_matrix(:,11,:)+par_matrix(:,25,:); %Activation module
    p2 = par_matrix(:,14,:)+par_matrix(:,26,:); % IkB module


    binscatter(p1(:),p2(:),100)

    hold on
    scatter(p1(:),p2(:),1,'filled')
    xlim([-20 10]);ylim([-20 10])

    hold on
    plot([-20 10],[-20,10],'--k','LineWidth',2)
    xlabel('Exit rate of free nuclear NFkB')
    ylabel('Exit rate of nuclear IkB-NFkB complex');

    subplot(2,2,3);

    p1 = par_matrix(:,18,:)+par_matrix(:,25,:); %Activation module
    p2 = par_matrix(:,16,:)+par_matrix(:,26,:); % IkB module


    binscatter(p1(:),p2(:),100)

    xlim([-20 10]);ylim([-20 10])

    hold on
    plot([-20 10],[-20,10],'--k','LineWidth',2)
    xlabel('IKK mediated degradation of bound IkB')
    ylabel('Basal degradation of bound IkB');

    subplot(2,2,4);

    p1 = par_matrix(:,17,:)+par_matrix(:,25,:); %Activation module
    p2 = par_matrix(:,15,:)+par_matrix(:,26,:); % IkB module


    binscatter(p1(:),p2(:),100)

    xlim([-20 10]);ylim([-20 10])

    hold on
    plot([-20 10],[-20,10],'--k','LineWidth',2)
    xlabel('IKK mediated degradation of free IkB')
    ylabel('Basal degradation of free IkB');
    %savefig(['param_covariation'],'png','-c0','-crop')
end