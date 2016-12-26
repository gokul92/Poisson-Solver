# Poisson-Solver

3-Dimensional serial Poisson solver has been developed using spectral methods (Fast-Fourier-Transforms) for spatial discretziation.
The solver has also been validated using standard test problems. The solution files are also present in the repository.
Separate subroutines/functions have been written to calculate the derivatives along x, y and z directions.
Derivatives are calculated in Fourier-Space as it is computationally inexpensive. FFTW library has been used to calculate the 
Fourier Transforms. Thus to run the code, it is necessary to have the library installed. 
Check the website http://www.fftw.org/fftw2_doc/fftw_6.html for more information on how to install the library.
