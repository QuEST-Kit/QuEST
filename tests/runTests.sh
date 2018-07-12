# grab the makefile
cp ../examples/makefile makefile

printf "============================\n"
printf "   COMPILING UNIT TESTS     \n"
printf "============================\n\n"

# compile QuEST and unit tests
make clean --silent
make EXE=runTests SOURCES=runTests COMPILER=gcc COMPILER_TYEP=GNU MULTITHREADED=0 DISTRIBUTED=0 GPUACCELERATED=0 --silent

# exit if compilation fails
if [ $? -eq 0 ]
then
	printf "\nCOMPILATION SUCCEEDED\n\n"
else
	printf "\nCOMPILATION FAILED\n\n"
	make clean --silent
	exit 1
fi

printf "============================\n"
printf "     RUNNING UNIT TESTS     \n"
printf "============================\n\n"

# run the unit tests
./runTests

# report test success
exitcode=$?
if [ $exitcode -eq 0 ]
then
	printf "\nTESTING PASSED\n\n"
else
	printf "\TESTOMG FAILED\n\n"
fi

# clean up and exit
make clean EXE=runTests --silent
rm makefile
exit $exitcode