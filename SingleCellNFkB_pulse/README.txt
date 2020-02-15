The data for this computational experiment was taken from Zhang et al 2017, Cell Systems.
Single-cell NFkB trajectories in response to 5 ng/ml 5 min. pulse TNF

The likelihood model assumes a standard deviation that is 10% of the data values (see config files)

(1) config_lasso, run_lasso and config_wo_lasso, run_wo_lasso are used for PT and PTLasso respectively. 
run_lasso initializes the PTLasso chains for all the NFkB trajectories with a previous short run that was fit to one of the trajectories.
run_wo_lasso randomly initializes the PT chains

(2) config_lasso_continue, run_lasso_continuechain and config_wo_lasso_continue, run_wo_lasso_continuechain are used for fixed start initializations from the data generated in (1)
