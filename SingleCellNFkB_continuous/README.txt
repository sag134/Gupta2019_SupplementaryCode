The data for this computational experiment was taken from Zhang et al 2017, Cell Systems.
Single-cell NFkB trajectory in response to 0.1 ng/ml continuous TNF stimulation.

run_lasso.m runs PTLasso from a randomized start point.
Here we use run_lasso.m using a large value of b (b=10, less constraint on the model so easier to fit) to get a short PT chain (3e5 samples)
The file obtained from run_lasso.m is used in run_lasso_fixedstart.m to initialize PTLasso with a smaller value of b (b=2, bigger constraint on the model so more difficult to fit) to get longer chains (9e5 samples).

call_screen_lasso.sh is used to start multiple screen sessions (each on a separate core), each running matlab in a headless mode with a different hyperparameter configuration.