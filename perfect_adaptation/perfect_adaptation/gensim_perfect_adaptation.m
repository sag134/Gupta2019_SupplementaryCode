function [err, sp, obsv] = gensim_perfect_adaptation(n,params)
%%
t = linspace(0,5,n)';
% run simulation
sp = [0,1,0];
obsv = [];
for i = 1:2
    sp(1) = i;
    [err, ~, sp, obsv_tmp] = perfect_adaptation( t, sp, params, 1 );
    sp = sp(end,:);
    if isempty(obsv)
        obsv = obsv_tmp;
    else
        obsv = [obsv;obsv_tmp(2:end)];
    end
end
if (err)
    sp = [];
    obsv = 1e29*ones(1,10);
    return;
end
return;

