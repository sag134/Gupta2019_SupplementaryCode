%This function calls an external function tight_subplot that can be
%obtained here: https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w

function [f,test_pt_no] = plot_post_separated(cfg,params_chain,energy_chain,param_indices,nbins,true_param,lim,LB,xlimits,frills,rl_index,adjust_y_axis)
%% function to plot posteriors from ptempest params_chain array for PTLasso examples
%% INPUTS:
% ptempest variables: cfg, params_chain, energy_chain
% param_indices: indices for plotting
% nbins: histogram bin number
% true_param: parameter values used to generate the exact model output
% lim: binary variable, if lim is 1 set the x axis limits to the limits
% from the cfg file
% xlimits: x axis limits for if lim is 0
% frills: Plot laplace prior boundaries if frills is 1
% rl_index: index of parameters used in the true data
% adjust_y_axis: Set the y axis the same for all plots. 
%% OUTPUTS:
% f: plot of the histograms, test_pt_no: number of points used per
% histogram

f = figure;
ub = find(energy_chain(1,:)==0,1)-1; % find how much of the energy chain is complete
if isempty(ub)
    ub = cfg.nswaps;
end
if isempty(param_indices)
    param_indices = 1:cfg.nparams;
end
if nbins == 0
    nbins = 1000;
end
c = 1;
l = @(x,b,mu)(1/(2*b))*exp(-abs(x-mu)/b); %laplace equation
M = numSubplots(length(param_indices));
m1 = M(1); m2 = M(2);
ha = tight_subplot(m1,m2,[.01 .01]);%,[.1 .01],[.01 .01])
test_pt_no = zeros(1,length(param_indices));
for i = param_indices
   % tight_subplot(m2,m1,c);
   axes(ha(c));
   c = c+1;
    %subplot(1,length(param_indices),c);c = c+1;
    p = params_chain(1,i,max(LB,cfg.adapt_last):ub); % get parameter chain
    hst = histogram(p(:),nbins,'EdgeColor','none','FaceAlpha',0.8);
    test_pt_no(i) = length(hst.Data);
    mxhst = max(hst.Values); %max histogram frequency
    hold on
    if isfield(cfg.param_defs{i},'min') && isfield(cfg.param_defs{i},'max') && lim == 1
        xlim([cfg.param_defs{i}.min cfg.param_defs{i}.max])
    elseif ~isempty(xlimits)
        xlim(xlimits)
    end
    ylim([0 mxhst]) % set y axis limits to max histogram frequency
    yl = ylim;
    % Plot true parameter values
    if ismember(i,rl_index)
        hold on
        plot([true_param(i),true_param(i)],yl,'--r','LineWidth',1)
        set(gca,'FontSize',15)
    end
    % Plot laplace prior boundaries
    if strcmp(cfg.param_defs{i}.prior,'laplace') || strcmp(cfg.param_defs{i}.prior,'boundedlaplace')
        b = cfg.param_defs{i}.b;
        mu = cfg.param_defs{i}.mu;           
        if isempty(xlimits)
            xl = xlim;
            xl(1) = min(xl(1),mu);
            xlim(xl);
        else
            xlim(xlimits)
            xl = xlim;
        end
        mx  = (1/(2*b));
        ll = mxhst*l(xl(1):0.0005:xl(2),b,mu)/(mx);
        hold on;
        %plot laplace equation
        if frills==1
            plot(xl(1):0.0005:xl(2),ll,'m','LineWidth',2);
        end
    end
end
%Adjust y axis if necessary
mx = 0;
mn = 0;
c = 1;
for i = 1:length(param_indices)
    %subplot(m2,m1,c);
    axes(ha(c));
    c=c+1;
    yl = ylim;
    mx = max(mx,yl(2));
    mn = min(mn,yl(1));
end
c=1;
for i = 1:length(param_indices)
    %subplot(m2,m1,c);
    axes(ha(c))
    c=c+1;
    if adjust_y_axis ==1
        ylim([mn,mx])
    end
    set(gca,'XTickLabels',[],'YTickLabels',[])
end
end