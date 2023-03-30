# SoMSP
This project is written in modern fortran and aimed at calculating precise values of scattering and extinction cross-sections and the scattering matrix for a single multilayered spheroidal particle.
It may be configured to use both double and native quadrouple precision.

The scattering equations are solved using the T-matrix method with expansion over the spheroidal functional basis. To calculate the spheroidal functions for complex parameter values we use
fortran subroutines created by Mr. Arnie Lee Van Buren (https://github.com/MathieuandSpheroidalWaveFunctions/).

The project is still in its testing phase. Use at your own peril, the feedback and bug reports are very welcome.

## Usage
### Building process
As a cmake project the code can be build in an empty directory with command 
```bash
cmake path_to_source_dir -DCMAKE_BUILD_TYPE={Release|Debug} [-DUSE_DOUBLE_PRECISION=True] [-DLOG={time;tmatrix;amplitude}]
```
then compile with
```bash
make main
```
and run
```bash
./main --input path_to_input_file [--scat-matr path_to_scattering_matrix_output]
```

### Building options
1. The spheroidal functions are calculated with quadrouple precision, but the general code for T-matrix solution can be configured to use double precision with -DUSE_DOUBLE_PRECISION=True. Otherwise 
quadrouple precision is used by default.
2. Option -DLOG= can be used to obtain the scattering.log file containing information of the last calculation. This option must be list of strings, separated
by symbol ';'. If it is empty, only the general information is printed. Additionally, tmatrices, amplitude matrices and times of execution of meaningful 
program blocks can be printed. For example, if times and tmatrices are required the build command would be:
```bash
cmake path_to_source_dir -DCMAKE_BUILD_TYPE=Release -DLOG=time;tmatrix
```
### Running options
1. The program requires either as input file specified with --input or a file input.txt in the current directory.
2. The output file for the scattering matrix can be specified with --scat-matr. If not specified, the matrix will be saved to scattering_matrix.txt.

### Input file format
It is a simple text file where each line contains a certain scatterer or light parameter. It must contain enough numbers, after which the rest of the string is ignored.
All values that depend on layer number are written from outside to the core.
1. f = -1|1 - indicated the type of the spheroid: -1 for oblate, 1 for prolate
2. n > 0 - number of layers
3. xv... - list of n real numbers - xv of the respective layer   
4. ab... - list of n real numbers - ab = a/b of the respective layer
5. ri... - list of n+1 complex numbers - refractive indices of the respective areas first number is the refractive index outside the particle, the rest inside its layers.
6. lamba - real number, the wavelength
7. alpha - the angle between the light propagation direction and the scatterer's symmetry axis in degrees.
8. lnum - the number of elements of expansion over the spheroidal functions. if it equals 0, the empirical value of 2*PI*a/lambda + 8 is used (here a is the bigger semiaxis of the spheroid)
9. spherical_lnum - depending on the calculation model sometimes our 'spheroidal' tmatrix is converted to a spherical basis. this parameter sets the size of that
'spherical' tmatrix. if it equals 0 it is set to be equal to lnum
10. minm - minimum value of m in the expansion. For a full calculation it should be 0.
11. maxm - maximum possible value of m in the expansion. The calculation may never reach this value of m as it stops when the next item in the sum becoes sufficiently
small in comparison with the already obtained result.
12. calculation model - string, expressing the method of calculation. it can have 1 of 5 values
- compare_uv_pq - calculate solution with UV potentials for m=0 and with PQ potentials. they must be equal and it's a good test that the general
spheroidal method works properly
- compare_te - calculate TM and TE mode using UV potentials and also TE from TM using the transition to spherical basis and back, as there is only one mode in
the spherical basis. Can be used to check the validity of the transition to the spherical basis
- all_uv - calculate for all m with UV potentials
- uv_pq - calculate for m > 0 with UV and additionally PQ replacing UV for m=0.
- uv_pq_te_from_tm - same as uv_pq, but the TE mode for m>0 in UV is obtained from TM

Generally, we consider uv_pq_te_from_tm to be the most precise and fastest way

Next 2 lines specify borders for scattering matrix calculation. All angles are in the coordinate system of the scatterer.

13. ntheta theta0 theta1 - calculate for ntheta uniformly distributed points in [theta0; theta1] (with ntheta-1 interval between them). 
theta0 and theta1 are in degrees and should be in [0; 360] with theta0 <= theta1
14. nphi phi0 phi1 same as theta, but phi0, phi1 should be in [0; 180]

Scattering matrix is calculated if both ntheta and nphi > 0.

The project contains as example file input.txt

### Building a library
The project has a library target:
```bash
make spheroidal_scattering
```
