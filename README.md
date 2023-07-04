corr_pc
-------

This is publicly available software for calculating correlation functions in galaxy survey data or
cosmological simulation snapshots (periodic boxes).  It has been used in practice for [Singh
(2021)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.508.1632S/abstract) and other works.


Installation
------------

After cloning the repository and entering its main directory, type `make` to compile it.  This
software requires the `mpic++` wrapper script that enables compiling MPI programs
([OpenMPI](https://www.open-mpi.org/)).

Usage
-----

`corr_pc` requires that the input data be in a particular format.  The notebooks in this repository
with names that start with `Gen_inp` provide documented examples of getting the input data into the
necessary format, including how to specify various parameters needed by the code to define the
correlation function calculations.  The notebooks require use of `astropy`.

The notebooks with names that start with `Process PC data` then provide examples of how to process
the input data and access the outputs.

License
-------

This software comes with no gurarantees/warranties and is solely provided in the hope that it may be
useful. If you use this code, please cite the [Singh
(2021)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.508.1632S/abstract) paper. You are free copy and
modify the code as desired as long a you provide appropriate citation to the original paper.



