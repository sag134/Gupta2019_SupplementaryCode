function f = plot_sim(cfg,params_chain,energy_chain,sim_handle,step,true_param,t,lg)
%%
    f = figure;
    ub = find(energy_chain(1,:)==0,1)-1; % find end point of energy chain
    if isempty(ub)
        ub = cfg.nswaps;
    end
    if step == 0
        step = 1000;
    end
    for i = cfg.adapt_last:step:ub
        i
        params = params_chain(1,:,i); % get parameter set
        [~,~, obsv] = sim_handle([],[], params);
        if isempty(t)
            t = linspace(1,cfg.data{1}.time(end),length(obsv));
        end
        if lg==1
            p = plot(t,log10(obsv),'b','LineWidth',2.5);
        else
            p = plot(t,obsv,'b','LineWidth',2.5);
        end
        p.Color(4) = 0.01;
        hold on
    end
    [~,~, obsv] = sim_handle([],[],true_param);
    if isempty(t)
        t = linspace(1,cfg.data{1}.time(end),length(obsv));
    end
    if lg==1
        plot(t,log10(obsv),'r','LineWidth',2.5);
     %  errorbar(cfg.data{1}.time,log10(cfg.data{1}.mean),log10(cfg.data{1}.stdev),'.r','LineWidth',2.5)
    else
        plot(t,obsv,'r','LineWidth',2.5);
        hold on
        errorbar(cfg.data{1}.time,cfg.data{1}.mean,cfg.data{1}.stdev,'.k','LineWidth',2.5)
        hold on
        scatter(cfg.data{1}.time,cfg.data{1}.mean,50,'k','filled');
    end
    hold on
    set(gca,'FontSize',15)
end