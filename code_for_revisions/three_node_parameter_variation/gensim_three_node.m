function [err, sp, obsv] = gensim_three_node(n,params)
t = linspace(0,50,n)';
% run simulation
[err, ~, sp, obsv] = three_node( t, [], params, 1 );
obsv = obsv(:,1);
if (err)
    sp = 1e29*ones(size(t));
    obsv = 1e29*ones(size(t));
    return;
end
return;

