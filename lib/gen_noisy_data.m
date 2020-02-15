function [smooth_data,noisy_data,mean_data,std_data,noise,f] = gen_noisy_data(true_param,simulate_function,noise_level,repeats,plt)
    %simulate data using true parameter values
    [~,~,smooth_data] = simulate_function(true_param);
    noisy_data = zeros(length(smooth_data),repeats);
    noise = zeros(length(smooth_data),repeats);
    %for each repeat
    for i = 1:repeats
        %the standard deviation of the noise is given by:
        noise(:,i) = (noise_level*smooth_data).*randn(length(smooth_data),1);
        noisy_data(:,i) = smooth_data+noise(:,i);
    end
    mean_data = mean(noisy_data,2);
    std_data = std(noisy_data,[],2); 
    f = 0;
    tpts = 1:length(smooth_data);
    if plt==1
        f = figure;
        plot(tpts,noisy_data);
        hold on
        errorbar(tpts,mean_data,std_data,'k','LineWidth',2)
        hold on
        plot(tpts,smooth_data,'--r','LineWidth',2); 
    end
end