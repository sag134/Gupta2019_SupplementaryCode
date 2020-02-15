function [err, sp, obsv] = simulate_three_node(t,init,params)
%%
% we are using the default initialization so init is left blank
t = linspace(0,50,8)';
% run simulation
[err, ~, sp, obsv] = three_node( t, [], params, 1 );
obsv = obsv(:,1);
if (err)
    sp = 1e29*ones(size(t));
    obsv = 1e29*ones(size(t));
    return;
end
return;

