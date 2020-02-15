This is the supplementary code for Gupta et al., "Parallel Tempering with Lasso for Model Reduction in Systems Biology"
For any updates to the code and/or documentation, please check here: https://github.com/RuleWorld/SupplementalMaterials/tree/master/Gupta2019

./<example_name> contains code to run each of the examples. key files in each of the example folders are:
    Model specification files: <example_name>.bngl, <example_name>.m, <example_name>_cvode.c
        the .bngl file is the model specification in bngl. 
        writeMfile and writeMexfile commands in the bngl file are used to get the .m file and the _cvode.c file
        the .c file is compiled to get a mex file.
        The simulate function calls the .m function which in turn calls the mex file for fast integration using CVODE solvers.

    <example_name>.mat 
        noisy data for fitting in a struct that can be accessed by the config file
    
    <example_name>_rawdata.mat
        contains the noise vectors used to get the noisy trajectories that are averaged for fitting

    simulate_<example_name>.m
        function that generates model output for fitting. Simulation has fixed temporal resolution. 
        Used to get 'observed data'

    gensim_<example_name>.m
        function that generates model output for fitting, but temporal resolution can be changed for visualizing results. 
        Used to get 'true data'

    config_lasso.m
        configuration file for parallel tempering with lasso using a bounded laplace prior

    config_wo_lasso.m
        configuration file for parallel tempering with uniform priors

    run_lasso.m
        function to run repeats of PTLasso for different hyperparameters

    run_wo_lasso.m
        function to run repeats of PT

    OPTIONAL FILES:Some additional optional files for easy scaling up:
            call_screen_lasso.sh
            call_screen_wolasso.sh
            run_matlab_with_lasso.sh
            run_matlab_wo_lasso.sh

./Figures
    contains code to generate each of the figures in the paper.
    All the code assumes that the data files by running each example are located
    in the corresponding example folder. 

./convergence_diagnostics
    code to check parameter convergence with PSRF or MPSRF, compute average step acceptance rates, 
    average swap acceptance rates and convergence of the energy chains for combining PT or PTLasso chains in the NFkB example
    All the PSRF and MPSRF calculations are done with third-party code available here: https://research.cs.aalto.fi/pml/software/mcmcdiag/

./lib

Descriptions of any additional example-specific files are provied in the relevant folders

code_for_revisions contains:
     (1) All the code for Figure S3.
     (2) reruns of the fits for NFkB responses to continuous TNF stimulation and calculation of the corresponding convergence metrics.
     (3) Any new figure code, or minor alterations to older figure code.
     (if code for a figure file is located in both FigureCode and code_for_revisions/revised_figure_code, the paper figure will correspond to the one in the latter.)


    