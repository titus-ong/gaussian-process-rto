====================
gaussian-process-rto
====================

Steady state real time optimisation using Gaussian processes to model systems

This repository contains MATLAB files for performing real time optimisation using Gaussian processes. The Gaussian process model is created using MATLAB's fitrgp function, and optimised using MATLAB's fmincon function. 

A simple optimisation problem, the Williams Otto reactor, is optimised in `example_WO.m </example_WO.m>`_. The model is shown to be applicable to more complicated flowsheets, exemplified in `example_hysys_fastrun.m </example_hysys_fastrun.m>`_ to optimise the HYSYS file `fast run.hsc </fast run.hsc>`_.

-----
To-Do
-----

MATLAB:

- Change TR radius with reduced gradient
- Save plots for excitation/no excitation, forgetting factor/no forgetting factor
- Expand to three variables (solvent flowrate, reboiler duty, inlet gas flowrate)

HYSYS:

- Level control storage drums