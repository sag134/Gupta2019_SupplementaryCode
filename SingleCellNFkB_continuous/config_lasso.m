function [cfg] = config_lasso(jobname,mu,b,trajectory_index)
% Initialize default configuration options 
cfg = init_config_defaults();

%% Parameter definitions (see BNGL model for parameter descriptions) [REQUIRED]
% Available prior distributions:
%   'point'           args: 'value'
%   'boundedlaplace','mu',MU,'b',B         args: 'min', 'max'
%   'normal',         args: 'min', 'sigma'
%   'lognormal',      args: 'min', 'sigma'
%   'boundedlaplace','mu',MU,'b',B,        args: 'min', 'max'
%   'boundeduniform', args: 'min', 'max', 'min', 'max'
%   'beta',           args: 'alpha', 'beta'
%   'boundedlaplace','mu',MU,'b',B,    args: 'min'


%LB and UB are lower and upper bounds for the unregularized
%reaction-specific parameters 
LB = -5;
UB = 10;
%LB1 and UB1 are lower and upper bounds for reguarization parameters
LB1 = mu-5;
UB1 = 6;
%Mu and b are the laplace parameters
B = b;
MU = mu;
%%
cfg.param_defs = { 
struct('name','Tot_NFKB' ,'prior','uniform', 'min', 1  , 'max',7 , 'units','uM' ), ...%NFKB protein
struct('name','IKK_N','prior','uniform', 'min', 1  , 'max',7,'units','uM' ), ...%IKK protein
struct('name','R',  'prior','uniform', 'min',1   , 'max',6 ,'units','uM' ), ...%Receptor protein
struct('name','k_b','prior','uniform','min',LB,'max',UB, 'units','uM'), ...% TNFR activated
struct('name','k_f','prior','uniform','min',LB,'max',UB, 'units','uM'), ...% basal deactivation of TNFR
struct('name','k_a', 'prior','uniform','min',LB,'max',UB, 'units','uM'), ...% IKK is activated
struct('name','k_4', 'prior','uniform','min',LB,'max',UB, 'units','uM'), ...% basal deactivation of IKK
struct('name','ka1a', 'prior','uniform','min',LB,'max',UB, 'units','uM'), ...% Complex formation
struct('name','kd1a', 'prior','uniform','min',LB,'max',UB, 'units','uM'), ...% Complex disassociation 
struct('name','ki1', 'prior','uniform','min',LB,'max',UB, 'units','uM'), ...% NFkB translocation into the nucleus
struct('name','ke1', 'prior','uniform','min',LB,'max',UB, 'units','uM'), ...% NFkB exit from nucleus
struct('name','ki2', 'prior','uniform','min',LB,'max',UB, 'units','uM'), ...% IKB entry into nucleus
struct('name','ke2', 'prior','uniform','min',LB,'max',UB, 'units','uM'), ...% IKB exit from nucleus
struct('name','ke2a','prior','uniform','min',LB,'max',UB, 'units','uM'), ...% IKB-NFKB complex exit from nucleus
struct('name','c4a','prior','uniform','min',LB,'max',UB, 'units','uM'), ...% Basal degradation of IKB
struct('name','c5a','prior','uniform','min',LB,'max',UB, 'units','uM'), ...% Basal degradation of IKB
struct('name','kt1a', 'prior','uniform','min',LB,'max',UB, 'units','uM'), ...% IKK mediated IkB degradation
struct('name','kt2a', 'prior','uniform','min',LB,'max',UB, 'units','uM'), ...% IKK mediated IkB degradation
struct('name','c1a', 'prior','uniform','min',LB,'max',UB, 'units','uM'), ...% IKB synthesis mediated by NFKB
struct('name','c_3', 'prior','uniform','min',LB,'max',UB, 'units','uM'), ...% A20 degradation
struct('name','c_1', 'prior','uniform','min',LB,'max',UB, 'units','uM'), ...% A20 synthesis mediated by NFKB
struct('name','a20_ikk_inact', 'prior','uniform','min',LB,'max',UB, 'units','uM'), ...%A20 mediated IKK inactivation
struct('name','a20_tnfr_inact', 'prior','uniform','min',LB,'max',UB, 'units','uM'), ...%A20 mediated TNFR inactivation
struct('name','A20_par', 'prior','boundedlaplace','mu',MU,'b',B,'min',LB1,'max',UB1, 'units','uM'), ...
struct('name','Activation_par', 'prior','boundedlaplace','mu',MU,'b',B,'min',LB1,'max',UB1, 'units','uM'), ...
struct('name','ikb_par', 'prior','boundedlaplace','mu',MU,'b',B,'min',LB1,'max',UB1, 'units','uM'), ...
};
%%
% initialize parameter distributions
cfg = init_parameter_defs( cfg.param_defs, cfg );
%% Observable definitions [REQUIRED]
%  Set display=0 to disable observable plot during visualization. Set
%  minplot/maxplot to see y-axis bounds in visualization scripts.
%  The 'norm' field should be the name of the observable that is used to
%  normalize the observable (e.g. YP is normalized by YT, the total quantity
%  of Y). Leave this field empty if normalization is not desired.
cfg.obsv_defs = {};
% initialize observable structs
cfg = init_observable_defs( cfg.obsv_defs, cfg );
%% Parallel temperating options
% Defaults are usually ok. Things you may want to change: jobname, nchains,
% parallel, nswaps, adapt_last, energy_init_max, relstep_init.
% See core/init_config_defaults.m for a complete list of config options.
cfg.jobname = jobname;               % job name, for file input/output
cfg.parallel = 0;                      % parallel? true/false
cfg.maxlabs  = 4;                      % maximinm number of labs for parallel computing
cfg.nchains  = 4;                      % number of chains
cfg.nswaps = 900000;                    % number of chain swaps (=number of saved samples!)
cfg.nsteps = 10;                       % number of stes between chain swaps
cfg.display_status_interval = 10;      % How often to display info
cfg.save_progress_interval = 50000;     % How often to save progress 
cfg.adapt_last = 50000;                 % last adaption step
cfg.energy_init_max = 1e12;%10000;             % maximinm allowed energy for initialization
cfg.beta_init = 0.6;                 % beta initialization parameter
cfg.relstep_init = 0.05;               % relstep initialization parameter
cfg.max_init_steps =5e10;% 5e6; 
%% Load experimental data [REQUIRED].
% The data file should import a cell array called "expt".
% Each cell array contains data for an experiment. It should define
%  the fields 'time', 'mean', 'stdev', 'nsamples' and 'weight'.
%  each array should have dimensions T x O, where T is the number of
%  time points and O is the number of observables, except for 'time'
%  which is is a column vector of length T. Set weight=0 if an observation
%  is missing or hidden (e.g. for future validation).  Missing data may be
%  indicated by NaN (not-a-number).
expt_data = load('data31_continuous.mat');   
expt = expt_data.expt;
% get number of experiments
cfg.nexpt = length(expt);
cfg.nexpt=1;
for d = 1 : cfg.nexpt
    % required fields for evaluating fit. The extra "0" added at the end is
    % supposed to be the value of IKKactive at the last time point, i.e.
    % IKK should adapt and activity should return to 0.
    % Here the requirement for IKK to adapt is a soft constraint (refer to
    % Gupta 2019 paper)
    cfg.data{d}.time = [expt{trajectory_index}.time,0];
    cfg.data{d}.mean = [expt{trajectory_index}.mean;0];
    cfg.data{d}.stdev = [0.1*expt{trajectory_index}.mean;0.1];
    cfg.data{d}.nsamples = [expt{trajectory_index}.nsamples;1];
    cfg.data{d}.weight   = [expt{trajectory_index}.weight;1];
end

%% load default functions
cfg = setup_default_functions( cfg );

%% setup custom function handles
% pass extra options to the protocol fcn here:
args = struct( ...
    'param_map',    cfg.param_map,    ...% useful for finding params by name
    'obsv_to_norm', cfg.obsv_to_norm, ...% required by "norm_obsv"
    'obsv_norm_by', cfg.obsv_norm_by  ...% required by "norm_obsv"
);
% experiment protocols here:
for d = 1 : cfg.nexpt
    % protocol specific parameters:
    % Experimental protocol function:
    %   prototype: [err,sim,obsv] = @(t,init,params)
    cfg.data{d}.protocol_fcn = @(t,init,params) simulate(t,init,params,args);

end
% Energy function [REQUIRED!]
%   prototype: [energy] = @(params)
%
cfg.energy_fcn = @(params) energy_gaussian(params, cfg);




