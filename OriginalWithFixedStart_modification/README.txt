to run ptempest from a fixed initial condition as in the case of the run files that point to OriginalWithFixedStart
replace this initialization file.

in the config file make the following changes:
(1) initialize cfg.initp to the desired initial parameter set
(2) set cfg.max_init_steps = -1 to prevent searching for initializations.