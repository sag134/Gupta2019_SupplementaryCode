function f = hyperparameter_tuning_plot(mean_gf,std_gf,mu_range,b_range)
%% 3D plot of mu,b,likelihood mean and standard deviation
    f = figure;
    for i=1:length(mu_range)
        plot3(mu_range(i)*ones(size(b_range)),log10(b_range),(mean_gf(i,:)),'LineWidth',2);%,std_gf,'LineWidth',2)
        hold on
        scatter3(mu_range(i)*ones(size(b_range)),log10(b_range),(mean_gf(i,:)),30,'filled','r');%,std_gf,'LineWidth',2)
        hold on
        for j = 1:length(b_range)
            j
            xV = [mu_range(i); mu_range(i)];
            yV = [log10(b_range(j)); log10(b_range(j))];
            zMin = mean_gf(i,j) + std_gf(i,j);
            zMax = mean_gf(i,j) - std_gf(i,j);
            zV = [zMin, zMax];
            % draw vertical error bar
            h=plot3(xV, yV, zV, '-k','LineWidth',1.5);hold on
            grid on
            xlabel('mu')
            ylabel('b')
            zlabel('-loglike')
            set(gca,'FontSize',22);
        end
    end
end

