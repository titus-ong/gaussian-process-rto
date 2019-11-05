====================
gaussian-process-rto
====================

Steady state real time optimisation using Gaussian processes to model systems

This repository contains MATLAB files for performing real time optimisation using Gaussian processes. The Gaussian process model is created using MATLAB's fitrgp function, and optimised using MATLAB's fmincon function. 

An example script `example_hysys.m </example_hysys.m>`_ is provided to show optimisation performed on the HYSYS file `testing.hsc </testing.hsc>`_.

-----
To-Do
-----

MATLAB:

- Forgetting factor for unconstrained problem
- Forgetting factor for expanding constraint (and shrinking)
- Normalise data
- Change TR radius with reduced gradient
- Put GP parameters into system files
- Change GP to model obj and constraints

HYSYS:

- Level control storage drums