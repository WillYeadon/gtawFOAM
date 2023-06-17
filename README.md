## gtawFOAM
OpenFOAM solver for Gas Tungsten Arc Welding

# Installation

I created this solver a while ago for OpenFOAM6. With some modifications it could be ported to later versions. Elsewise you'll need to install Openfoam6 and put this in the same directory as your custom solvers. Put all the code into a folder called 'gtawFOAM' in your custom solver directory and run wclean & wmake. The solver should then run with the command gtawFOAM. See the repo 'thesis' for the cases I used. You can read the full thesis [here](https://arxiv.org/abs/2205.07687). Elsewise, looking at other repos it appears there is one called [beamWeldFoam](https://github.com/tomflint22/beamWeldFoam) that may be of interest if you're looking at this solver.
