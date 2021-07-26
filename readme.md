# 2021 Renewal Equations

This folder contains some codes that are useful to reproduce the results in the paper
Scarabel F, Diekmann O, Vermiglio R (2021). Numerical bifurcation analysis of renewal equations via pseudospectral approximation, 
Journal of Computational and Applied Mathematics, 397:113611. https://doi.org/10.1016/j.cam.2021.113611
Please cite the reference when using this code.

The codes are used to perform the numerical bifurcation analysis of the examples in the paper.
The numerical bifurcation is performed using the Matlab package MatCont, available at: https://sourceforge.net/projects/matcont/

The interpolation is performed using polint.m contained in the Differentiation Matrix Suite
from Weideman, Reddy 2000: http://appliedmaths.sun.ac.za/~weideman/research/differ.html

For the MatCont continuation, each example consists of two files:
1) PS_example: matlab function containing the definition of the right-hand side of the ODE system obtained through pseudospectral discretization, in the format suitable for the Matcont continuation.
2) MC_example: script for the MatCont continuation of the system defined in "PS_example".

The codes are tested on MATLAB 2020b and MatCont version MatCont7p1.
