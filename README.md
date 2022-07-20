# commSelect
Simulates artificial selection of microbial communities.

Dependencies for this package are NumPy v1.18.0 and SciPy v1.4.1.  To run the mature step in parallel, you will also need the Python package mpi4py.futures.

To install, simply clone this repository in your directory of choice.  You can test the package by running HMSixPhenotypes.py. New files that use commSelect should be placed in the same directory as commSelect/.

 During the simulation, a random number generator `rng_main` is pickled at the beginning of each cycle. It is then used to generate `num_wells` random integers, each serving as the random number generator seed for each community.

To replay a certain cycle or community, load the `rng_main` for that cycle. For example, one can run `HMSixPhenotypes.py` first. The results for Cycle 0 and 1 would be saved. In `replay.py`, the random number generator and chosen Adults are loaded from Cycle 0. Cycle 1 is then replayed. The replayed results are the same as the saved results.

My modification:
* Implement the current random number interface (e.g., instead of using `np.random.seed` to set a global random stream, use `np.random.default_rng` to pass around a random stream)
* Assign a random number generator seed for each community so that the maturation can be replayed.
* The random number generator for reproduction is saved for each cycle so that any cycle can be replayed without going through the previous cycles.
* Data are saved through pickling instead of as a text file.
