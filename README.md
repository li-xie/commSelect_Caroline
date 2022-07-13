# commSelect
Simulates artificial selection of microbial communities.

Dependencies for this package are NumPy v1.18.0 and SciPy v1.4.1.  To run the mature step in parallel, you will also need the Python package mpi4py.futures.

To install, simply clone this repository in your directory of choice.  You can test the package by running HMSixPhenotypes.py. New files that use commSelect should be placed in the same directory as commSelect/.

My modification:
* Implement the current random number interface (e.g., instead of using `np.random.seed` to set a global random stream, use `np.random.default_rng` to pass around a random stream)
* Assign a random number seed for each community so that the maturation can be replayed.
* The random number generator for reproduction is saved for each cycle so that any cycle can be replayed without going through the previous cycles. 

