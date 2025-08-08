module load profile/candidate
module load nvhpc/25.3
cp Makefile_leonardo Makefile
rm *.mod
rm mhit36
make clean
make

rm -rf output
mkdir output
rm *.o