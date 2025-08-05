NVARCH=Linux_x86_64; export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/25.7/compilers/man; export MANPATH
PATH=$NVCOMPILERS/$NVARCH/25.7/compilers/bin:$PATH; export PATH
export PATH=$NVCOMPILERS/$NVARCH/25.7/comm_libs/mpi/bin:$PATH
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/25.7/comm_libs/mpi/man
cp Makefile_local Makefile
rm *.mod
rm nemesi36
#nvfortran -fast  -gpu=managed -acc -cudalib module.f90  poissonfast.f90 main.f90 -o mhit36 -L/usr/local/cuda/lib64 -lcufft
make

rm -rf output
mkdir output
rm *.o
rm *.dat

#run the code
./nemesi36
