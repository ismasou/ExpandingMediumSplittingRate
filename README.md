# Splitting Rate Solver in Expanding medium

This is a simple solver for QCD in-medium splitting rate in evolving media.
The Makefile contains few options to compile Static and Expanding medium solvers.
To get a feel of how to run the code use the bash scripts `CH_G_Check.sh` which generates the splitting rate in a static medium.
It will also generate a plot to check against [1] results in the `CH_G_Comparison` folder.

For an expanding medium, use the `./Expanding.sh` script.
This will compile the code and generate the splitting rates and plots for the expanding medium in `ExpandingPlot` folder.

## Usage

There are several executables that can be compiled:
    - `Static`: for a static medium
    - `Expanding`: for an expanding medium
    - `Opacity`: for opacity expansion in static medium
    - `ExpandingOp`: for opacity expansion in expanding medium

All the executables take the following arguments:
```bash
    Usage
    -T: Temperature
    -P: Parent momentum
    -z: Momentum fraction
    -tmax: Maximum time
    -tmaxfm: Maximum time in fm
    -t0: Starting time of Bjorken expansion
    -a: Exponent
    -t1: The time the hard parton enters the medium
```

## Dependencies
This code uses the following libraries:
    - `gsl`: gnu scientific library.
    - `fmt`: formatting library for output.
    - `openmp`: for parallelization.
    - `boost/math/quadrature/gauss_kronrod.hpp`: for numerical integration.


## References
[1] S. Caron-Huot and C. Gale, Phys. Rev. C 82, 064902 (2010), arXiv:1006.2379 [hep-ph].
