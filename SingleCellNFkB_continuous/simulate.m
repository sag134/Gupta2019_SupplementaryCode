function [err, sp, obsv] = simulate( t, init, parameters,cfg)
num_par = 26;
params = parameters(1:num_par);
params = [params,0];
tequil = (0:1e5:1e8)';
tstim = (0:5*60:(48-1)*5*60)';
dose = 0.1;
t0 = 4;
species_init = [];
t = [tstim;1];
%% add hard constraints on parameter covariation
% (1) params(14) = exit rate of nuclear bound NFkB-IkB. This should be
% greater than exit rate of free IkB (13) or (2) free NFkB (11)
% (3) Basal degradation of bound IkB (c5a,16) should be less than IKK mediated
% degradation rate of bound IkB (kt2a,18)
% (4) Basal degradation rate of free IkB (c4a, 15) should be less than IKK
% mediated degradation rate of free IkB (kt1a,17)
if (params(14)+params(26))<(params(13)+params(26)) || (params(14)+params(26))<(params(11)+params(25)) || (params(16)+params(26))> (params(18)+params(25)) || (params(15)+params(26))>(params(17)+params(25))
    err = 1;
    sp = 1e29*ones(size(t));
    obsv = 1e29*ones(size(t));
    return
end
%% Equilibration
[err, ~, eq_species_out,oeq] = nfkb_lasso(tequil,species_init,params,1);
if err
    sp = 1e29*ones(size(t));
    obsv = 1e29*ones(size(t));
    return
end
%Check absolute tolerance for equilibration
if max((oeq(end,:)-oeq(end-5,:))) >=1e-4
    sp = 1e29*ones(size(t));
    obsv = 1e29*ones(size(t));
    return
end
%% Check equilibrium fraction of NFkB in the nucleus. Set a wide allowable range of 5-50%.  
NuclearAbundance = (oeq(end,1)+oeq(end,6))/10^params(1);
if NuclearAbundance<0.05 || NuclearAbundance > 0.50 
    err = 1;
end
if err
    sp = 1e29*ones(size(t));
    obsv = 1e29*ones(size(t));
    return
end
%% TNF stimulus    
species_perturb = eq_species_out(end,:); species_perturb(end) = dose;        
[err, ~, species_out, o1] = nfkb_lasso(tstim,species_perturb,params,1);
if err
    sp = 1e29*ones(size(t));
    obsv = 1e29*ones(size(t));
    return
end    
%% Check numbers of active ikk and tnfr
ikk = o1(:,8);
tnfr = o1(:,9);
if max(ikk) < 1 || max(tnfr) < 1 
    err = 1;
end
if err
    sp = 1e29*ones(size(t));
    obsv = 1e29*ones(size(t));
    return
end    
obsv = [oeq(end-t0+1:end,1)+oeq(end-t0+1:end,6);o1(1:end-t0,1)+o1(1:end-t0,6)]/(oeq(end,1)+oeq(end,6));
obsv =[ obsv;ikk(end)];
sp = [];
return    
end  
