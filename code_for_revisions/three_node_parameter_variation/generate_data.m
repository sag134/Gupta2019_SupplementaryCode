addpath('/shared2/LabUserFiles/Sanjana_Gupta/LassoManuscript/lib')

for fr_error = 0.05
    f = figure('Units','normalized','OuterPosition',[0,0,1,0.5]);
    num_repeats=10;
    parameters = [ -1, -10, 0, -10, -10, -10 ]; 
    [err, sp, obsv] = simulate_three_node([],[],parameters);
    smooth_data=obsv;
    subplot(3,6,1:3);
    plot(obsv,'r','LineWidth',2);
    data=zeros(size(obsv,1),num_repeats);
    params=zeros(num_repeats,length(parameters));
    hold on
    for i=1:num_repeats
        params(i,:) = parameters;
        params(i,parameters~=-10) = parameters(parameters~=-10) + fr_error*randn(1,length(parameters(parameters~=-10)));
        [err, sp, obsv] = simulate_three_node([],[],params(i,:));
        data(:,i) = obsv;
        plot(obsv);
    end

    subplot(3,6,4:6);
    errorbar(1:8,mean(data,2),std(data,[],2));
    for i = 1:6
        subplot(3,6,i+6);
        histogram(params(:,i));
        yl = ylim;
        hold on
        plot([parameters(i),parameters(i)],yl,'r');
        xlim([-12 3])
        subplot(3,6,i+12);
        histogram(params(:,i));
        yl = ylim;
        hold on
        plot([parameters(i),parameters(i)],yl,'r');

    end
    mean_data = mean(data,2);
    std_data = std(data,[],2);
    expt = get_expt_struct(mean_data,std_data);

 %   save(['three_node_parameter_variability',num2str(fr_error),'.mat'],'expt')
 %   save(['full_three_node_parameter_variability',num2str(fr_error),'.mat'])
 %   savefig(['three_node_parameter_variability',num2str(fr_error*100)],'png','-c0','-crop')

end