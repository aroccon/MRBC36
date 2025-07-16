## NEMESI36
2D code for PFM simulations (WIP)


## Check list of features implemented (or to be implemented)

- Poisson solver (Matlab) ✅
- Poisson solver validation (Matlab) ✅
- Poisson solver convergence (Matlab) ✅
- Poisson solver (Fortran: FFT forth and back) ✅
- Poisson solver (Fortran) ✅
- Added support for stretched nodes along z (Matlab only, easy to move in Fortan) ✅
- Introduced phi ✅ 
- Introduced temperature ✅ 
- Introduced temporal loop ✅ 
- Navier-Stokes solution ✅ 
- Validation (phi) ✅ 
- Validation (temp) ✅ 
- Validation (NS): Channel flow (laminar + moving top wall) ✅
- RB setup ✅ 
- GPU offloading of entire code  ✅
- TDMA optimization (x 10 speed-up)  ✅

## RBC (Ra=1e8)

![Test](doc/rbc2.png)


## Grid points (staggered)

![Test](doc/grid.png)


