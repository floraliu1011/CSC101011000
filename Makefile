SHELL := /bin/bash
FORWARD_PATH=test_forward.txt#fill this in#
REVERSE_PATH=test_reverse.txt#fill this in#
NUM_FILE=4

total=`wc -l ${FORWARD_PATH} | awk '{print $1}'`
if ((${NUM_FILE} - 1)) && (( ${total} % (4 * (${NUM_FILE} - 1)) )); then
				line=$(( 4 *(${total} / (4 * (${NUM_FILE} - 1))) ))
else
				line=$(( 4*(${total} / (4 * (${NUM_FILE})) ))
fi
FORWARD_NAME=`cut -d "." -f1 ${FORWARD_PATH}`
REVERSE_NAME=`cut -d "." -f1 ${REVERSE_PATH}`

submit: submission
	qsub submission

submission:
	cat > submission <<EOF
	#!/bin/bash
	#PBS -l nodes=$(( $NUM_FILE / 16 + 1)):ppn=8,walltime=6:00:00
	#PBS -o output.out
	#BPS -e output.err
	
	module load gnu-parallel/20140622
	module load intel/15.0.2
	module load python/3.5.1
	
	cd $PBS_O_WORKDIR
	EOF
	for i in seq ${NUM_FILE}
	do
		echo "(python parser.py Intermediate/${FORWARD_NAME}$i Intermediate/${REVERSE_NAME}$i Output/LogFile/log_file$i.txt)&" >> submission
		echo wait >> submission
	done

prepare:
	mkdir Data Intermediate Output Output/LogFile
	mv ${FORWARD_PATH} ${REVERSE_PATH} Data
	echo finished preparing

split: ${FORWARD_PATH} ${REVERSE_PATH}
	split -a $(( `echo $NUM_FILE | wc -m` - 1 )) -l ${line} -d ${FORWARD_PATH} ${FORWARD_NAME}
	split -a $(( `echo $NUM_FILE | wc -m` - 1 )) -l ${line} -d ${REVERSE_PATH} ${REVERSE_NAME}
	mv ${FORWARD_NAME}* Intermediate
	mv ${REVERSE_NAME}* Intermediate
	echo finished splitting
