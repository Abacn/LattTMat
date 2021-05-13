# Transfer matrix method for next-nearest-neighbor Ising models

## Background

The transfer matrix (TM) method calculate the partition function and correlation length of a semi-infinite stripe of $L \times \infty$, obtained from the first a few eigenvalues. It also gives the marginal probability distribution of every possible state on a layer $L$. Combining these one can also generate the equilibrium configuration of strip width L and arbitrary length N.

## Summary

Language: C++ (std=c++11)
Publications:
- Yi Hu and Patrick Charbonneau. "Resolving the two-dimensional axial next-nearest-neighbor Ising model using transfer matrices." Physical Review B, 103, 094441 (2021). http://doi.org/10.1103/PhysRevB.103.094441

## File descriptions

- TMat/ : transfer matrix computation core
- DumpConfig/ : dump equilibrium configurations
- tmat.make : makefile script for transfer matrix calculation
- dumpconfig.make : makefile script for configuration planting


## Compilation

The code requires [Eigen3](https://gitlab.com/libeigen/eigen.git)

The code also requires [Spectralib v0.9](https://github.com/yixuan/spectra/releases/tag/v0.9.0). Note that v1.x  the API has changed.

Install these libs in system, or set the header path in $INCLUDE environment variable. Run `make -f tmat.make` to compile for tmat and `make -f dumpconfig.make` for dumpconfig. Modify makefile scripts when necessary.


## Running

### Transfer matrix calculation

The binary "tmat" takes one parameter which is the input file. It computes for the eigenvalue, right and left eigenvectors output as "X_val.dat", "X_vec.dat", "X_vecl.dat". 

An example input file is illustrated below:

```
modeltype = B # model type. A: ANNNI; B: BNNNI and 3NN; D: DNNI
dim = 2       # dimensionality. 2 or 3
axialdir = 6  # indices of the type of boundary condition (see source)
neigs = 1     # number of eigenvalues been calculated
savevec = 1 1 # options (1 or 0) to save vectors or computing layer magnetization (or not)
outfile = 2DYBM  # prefix for output files
size = 12        # system width L
J = 0.5          # inverse temperature beta * J 
kappa = 0.8 0    # kappa1 and kappa2 values
h = 0            # external field
```

For 2D (dim = 2):
- modeltype support A, B, D for the ANNNI, BNNNI+3NN, and DNNI model, respectively.

- For the ANNNI model, axialdir supports
-- 0: y-double-layer
-- 1: $\perp$-TM, full $2^L \times 2^L$ transfer matrix
-- 2: $\parallel$-TM, full $4^L \times 4^L$ transfer matrix
-- 5: $\perp$-TM, reduced size
-- 6: $\parallel$-TM, reduced size
-- 9: $\perp$-TM, open boundary
-- 13: $\perp$-TM, reduced size, open boundary
-- 17, 18: $\perp$-TM, a pair of additional up-up or up-down spins in each side
-- 19, 20: $\perp$-TM, two pair of additional up-up or up-down spins in each side

- For the BNNNI/3NN model, axialdir supports
-- 3: $/$-TM, full $2^L \times 2^L$ transfer matrix
-- 6: $\parallel$-TM, reduced size
-- 7: $/$-TM, open boundary
-- 13: $\perp$-TM, reduced size

- For the DNNI model, axialdir supports
-- 1: $\perp$-TM, full $2^L \times 2^L$ transfer matrix
-- 2: $\parallel$-TM, full $4^L \times 4^L$ transfer matrix
-- 5: $\perp$-TM, reduced size
-- 6: $\parallel$-TM, reduced size
-- 17, 18: $\perp$-TM, a pair of additional up-up or up-down spins in each side

### Configuration planting

The binary "dumpconfig" takes three parameters, for example, "./dumpconfig input_ymannni.dat 1 100". The first parameter is the input file, the second is the number of configurations to generate, the third is the row of configurations.

It generates an equilibrium configuration, e.g.,

```
000000000000
000000000000
111111111111
111111111111
011100000000
000000000000
100011111111
111111111111
111111111111
000000000000
```

where 0s are up spins and 1s are down spins.



