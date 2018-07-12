# user settings for test compilation
COMPILER=gcc
COMPILER_TYPE=GNU
MULTITHREADED=1
DISTRIBUTED=0
GPUACCELERATED=0
GPU_COMPUTE_CAPABILITY=30

# grab the makefile
cp ../examples/makefile makefile

printf "============================\n"
printf "   COMPILING UNIT TESTS     \n"
printf "============================\n\n"

# compile QuEST and unit tests
make clean --silent
make EXE=runTests SOURCES=runTests COMPILER=$COMPILER COMPILER_TYPE=$COMPILER_TYPE MULTITHREADED=$MULTITHREADED DISTRIBUTED=$DISTRIBUTED GPUACCELERATED=$GPUACCELERATED GPU_COMPUTE_CAPABILITY=$GPU_COMPUTE_CAPABILITY --silent

# exit if compilation fails
if [ $? -eq 0 ]
then
	printf "\nCOMPILATION SUCCEEDED"
	printf "\n---------------------\n\n"
else
	printf "\nCOMPILATION FAILED"
	printf "\n---------------------\n\n"
	make clean --silent
	exit 1
fi

printf "============================\n"
printf "     RUNNING UNIT TESTS     \n"
printf "============================\n\n"

# run the unit tests
export OMP_NUM_THREADS=3
if [ $DISTRIBUTED == 0 ]
then
    ./runTests
else
    mpirun $MPI_HOSTS ./runTests
fi

# report test success
exitcode=$?
if [ $exitcode -eq 0 ]
then
	printf "\nTESTING PASSED"
	printf "\n--------------\n\n"
else
	printf "\nTESTING FAILED"
	printf "\n--------------\n\n"
fi

# clean up and exit
make clean EXE=runTests --silent
rm makefile
exit $exitcode
