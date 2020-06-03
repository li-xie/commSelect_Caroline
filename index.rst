.. BioCommSelect documentation master file, created by
   sphinx-quickstart on Thu May 14 14:21:23 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to BioCommSelect's documentation!
=========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

--------
Overview
--------

BioCommSelect is a Python package that simulates artificial selection of microbial communities. It was created for the Shou group.

Artificial selection is used to "optimize" a complex ecosystem by creating many replicates of that ecosystem and selecting the most "optimal" (best at performing some desired function) communities for reproduction after some period of time. This can counteract the process of natural selection moving a biological community away from its desired performance in exchange for improved fitness.

The details of an artificial selection protocol have a large impact on whether or not the protocol will be successful. Considerations include the number of communities per "cycle", the time they are given to mature, the method of picking the "best" communities, the method of propagating those communities to new wells, and more. This package is designed for biologists to "test out" different protocols before performing them in a wetlab setting.

--------------------
Defining a Community
--------------------

In bioCommSelect, a community is defined by a set of "cell types" (distinct species of organisms with unique traits that may mutate) and metabolites. A community's dynamics are defined by a system of differential equations describing (1) the rate of change of relevant metabolite concentrations, and (2) the birth rate per biomass of each cell strain. The user initializes the cell types and metabolites and writes a function containing the system of differential equations.  From there, bioCommSelect can run the "maturation" phase, where the community grows and mutates.  Mature communities can then be ranked based on their community function, and the best ones are pipetted into new wells to perform the maturation phase again.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
