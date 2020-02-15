function [err, sp, obsv] = simulate_linear_full_trajectory(d,t,init,params)
%%
t = linspace(0,800,4)';
% run simulation
sp = [d,0,0];
[err, ~, sp, obsv] = perfect_adaptation( t, sp, params, 1 );
if (err)
    sp = [];
    obsv = 1e29*ones(1,length(t));
    return;
end
return;

