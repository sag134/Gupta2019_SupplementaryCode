function [err, sp, obsv] = simulate_five_node(t,init,params)
%%
t = linspace(0,50,8)';
% run simulation
[err, ~, sp, obsv] = five_node( t, [], params, 1 );
obsv = obsv(:,1);
if (err)
    sp = 1e29*ones(size(t));
    obsv = 1e29*ones(size(t));
    return;
end
return;

