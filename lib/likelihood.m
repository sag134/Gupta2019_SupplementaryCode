function [energy] = likelihood( params, cfg )
%This is equivalent to the energy_gaussian_wopenalty function

energy = 0;%-logprior;
%%
% start integration timer
simtimer = tic;

% % equillibrate, if required
if isfield(cfg, 'equilibrate_fcn')
    [err,state,~] = cfg.equilibrate_fcn( params );
    if (err)
        energy = cfg.big_energy;
        return;
    end
    init = state(end,:);
else
    init = [];
end

% simulate experiments
for d = 1:cfg.nexpt

    % simulate experiment (default initial conditions)
    [err,~,obsv] = cfg.data{d}.protocol_fcn( cfg.data{d}.time, init, params );
    if (err)
        energy = cfg.big_energy;
        return;
    end

    % normalize obsv
  %  [obsv] = norm_obsv( obsv, params, cfg );

    % heuristic penalties (do this before transformating obsv)
    if isfield(cfg.data{d}, 'heuristic_penalty_fcn')
        penalty = cfg.data{d}.heuristic_penalty_fcn(obsv, params);
        if isinf(penalty)
            energy = cfg.big_energy;
            return;
        end
        energy = energy + penalty;
    end

    % if necessary, transform simulated trajectory for computing fitness
    if isfield(cfg, 'transform_sim_for_fit_fcn')
        obsv = cfg.transform_sim_for_fit_fcn(obsv,params);
    end

    % calculate log-likelihood as weighted sum of square errors
    loglike = nansum(nansum( -cfg.data{d}.weight .* (obsv - cfg.data{d}.mean).^2 ./ (2*(cfg.data{d}.stdev).^2) ));
    % subtract likelihood from energy
    energy = energy - loglike;

end

% penalize for slow integrations
dt = toc(simtimer);
energy = energy + cfg.timepenalty*dt^2;

% all done
return;
