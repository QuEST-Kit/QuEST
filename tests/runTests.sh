
# grab the makefile, use its current settings
cp ../makefile makefile

printf "============================\n"
printf "   COMPILING UNIT TESTS     \n"
printf "============================\n\n"

# compile QuEST and unit tests
make clean --silent
make EXE=runTests SOURCES=runTests QUEST_DIR=../QuEST --silent

# exit if compilation fails
if [ $? -eq 0 ]
then
	printf "\nCOMPILATION SUCCEEDED"
	printf "\n---------------------\n\n"
else
	printf "\nCOMPILATION FAILED"
	printf "\n---------------------\n\n"
	make clean EXE=runTests --silent
	exit 1
fi

printf "============================\n"
printf "     RUNNING UNIT TESTS     \n"
printf "============================\n\n"

# run the unit tests
export OMP_NUM_THREADS=3
distributed=$(make SUPPRESS_WARNING=1 getvalue-DISTRIBUTED SILENT=1 --silent)
if [ $distributed == 0 ]
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
