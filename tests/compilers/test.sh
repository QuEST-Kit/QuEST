# Tyson Jones, 10th July
# tyson.jones@ccc.ox.ac.uk
# 
# Automates testing of QuEST with different C/C++/NVCC/MPICC compiler versions
# Requires the following files which contain newline-separated lists of the names
# of SLURM modules to load:
# - gcc_versions.txt
# - icc_versions
# - cuda_versions
# - mpi_versions

sourcename="testcode"
execname="testexec"
outputname="compatibility.csv"

declare -a compiler_types=("GNU" "INTEL")
declare -a wrappers=("NONE" "MPICC" "NVCC")
declare -a languages=("c" "cpp")

# write table headers to file
echo "LANGUAGE,COMPILER,VERSION,WRAPPER,VERSION,COMPILED,RAN\n" > $outputname

# try every language
for language in "${languages[@]}"
do
	# clear file extension
	for reset in "${languages[@]}"
	do
		mv "$sourcename.$reset" "$sourcename" 2>/dev/null
	done
	
	# set file extension
	mv "$sourcename" "$sourcename.$language"

	# try every type of compiler
	for compiler_type in "${compiler_types[@]}"
	do
		# (determine compiler cmd)
		if [ $compiler_type == "CLANG" ]
		then
			compiler="clang"
		elif [ $compiler_type == "GNU" ]
		then
			compiler="gcc"
		elif [ $compiler_type == "INTEL" ]
		then
			compiler="icc"
		else
			printf "\nERROR 0\n$compiler_type\n"
			exit 125
		fi
	
		# get compiler versions
		readarray -t compiler_versions < "${compiler}_versions.txt"
	
		# try every type of wrapper
		for wrapper in "${wrappers[@]}"
		do
			# determine platform args from wrapper
			if [ $wrapper == "NONE" ]
			then
				distributed=0
				gpuaccelerated=0
				declare -a wrapper_versions=("")
			elif [ $wrapper == "MPICC" ]
			then
				distributed=1
				gpuaccelerated=0
				readarray -t wrapper_versions < "mpi_versions.txt"
			elif [ $wrapper == "NVCC" ]
			then
				distributed=0
				gpuaccelerated=1
				readarray -t wrapper_versions < "cuda_versions.txt"
			else
				printf "\nERROR 1\n$wrapper\n"
				exit 125
			fi
			
			# try every compiler version
			for compiler_version in "${compiler_versions[@]}"
			do
				# try every wrapper version
				for wrapper_version in "${wrapper_versions[@]}"
				do
					
					printf "$language: ($compiler_version) wrapped by ($wrapper_version)\n"
					
					# test loading compilers
					module purge
					result=$(module load $compiler_version $wrapper_version 2>&1 >/dev/null)
					
					# skip incompatible compilers
					if [ ! -z  "$result" ]
					then
						printf "\nINCOMPATIBLE\n\n"
						continue
					fi
					
					# actually load compilers
					module purge
					module load binutils $compiler_version $wrapper_version 2>&1 >/dev/null
					
					# test compilation
					make clean --silent
					make COMPILER=$compiler COMPILER_TYPE=$compiler_type DISTRIBUTED=$distributed GPUACCELERATED=$gpuaccelerated --silent
					
					# record compilation success
					if [ $? -eq 0 ]
					then
						printf "\nCOMPILATION SUCCEEDED\n\n"
						comp_success=1
					else
						printf "\nCOMPILATION FAILED\n\n"
						comp_success=0
					fi
					
					# test and record execution success
					if [ $comp_success -eq 1 ]
					then
						printf "Executing...\n"
						./testexec
						if [ $? -eq 0 ]
						then
							printf "\nEXECUTION SUCCEEDED\n\n"
							exec_success=1
						else
							printf "\nEXECUTION FAILED\n\n"
							exec_success=0
						fi
					else
						exec_success=0
					fi
					
					# replace NONE wrapper with emptry string (prettier)
					wrapper_str=$wrapper
					if [ $wrapper_str == "NONE" ]
					then
						wrapper_str=""
					fi
					
					# append data to file
					row="$language,$compiler_type,$compiler_version,$wrapper_str,$wrapper_version,$comp_success,$exec_success"
					echo $row >> $outputname
					
				done
			done
		done
	done
done

# clean up the final compilation
make clean
