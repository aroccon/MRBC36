NVARCH=Linux_x86_64; export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/24.3/compilers/man; export MANPATH
PATH=$NVCOMPILERS/$NVARCH/24.3/compilers/bin:$PATH; export PATH
export PATH=$NVCOMPILERS/$NVARCH/24.3/comm_libs/mpi/bin:$PATH
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/24.3/comm_libs/mpi/man
cp Makefile_local Makefile
rm *.mod
rm nemesi36
make

rm -rf output
mkdir output
rm *.o
rm *.dat

#run the code
nsys profile --trace=cuda,osrt,nvtx,openacc ./nemesi36
