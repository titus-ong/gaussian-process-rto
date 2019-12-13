====================
gaussian-process-rto
====================

Steady state real time optimisation using Gaussian processes as MA modifiers

This repository contains MATLAB files for performing real time optimisation using Gaussian processes. The Gaussian process model is created using MATLAB's fitrgp function, and optimised using MATLAB's fmincon function. 

A simple optimisation problem, the Williams Otto reactor, is optimised in `example_WO.m </example_WO.m>`_. The model is shown to be applicable to more complicated flowsheets, exemplified in the various example_hysys files to optimise HYSYS flowsheet models of the Imperial College Carbon Capture Pilot Plant (IC CCPP). Each example file is a different optimisation problem.

---
Use
---

1. Adjust system settings by editing the system files:

   - Tolerances
   - Constraints
   - Starting point
   - Trust region parameters

2. Set excitation variable (GP.excite) to true or false in example file
3. Run the example file to solve that particular optimisation problem. E.g. `example_hysys_fastrun.m </example_hysys_fastrun.m>`_ is a 2-input variable economic optimisation problem.

Operating plots of the GP model can also be viewed by running the relevant script found in `/System Files/Operating region plots </System Files/Operating region plots>`_.