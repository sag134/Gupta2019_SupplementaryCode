./perfect_adaptation contains files for running PT and PTLasso from
a randomized start point. However since it takes a while to find the minima,
more than half of the chains consisted of burn-in samples that correspond to bad fits

./perfect_adaptation_continue_chain contains files that use data from the ./perfect_adaptation folder
to initialize PTLasso and PT chains with a parameter set from the minima. This way the full chain is sampling the desired minima
instead of having to discard the a large fraction of the chain as burn-in.

This process can be repeated if necessary (if some chains are still not converging in the set number of samples), and
./perfect_adaptation_continue_chain_two has files to run a third set of chains intialized from data in ./perfect_adaptation_continue_chain

