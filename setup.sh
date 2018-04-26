#!/bin/bash
set -e #abort the process if encountered error
# check cmdline args
if test $# -ne 3
then
	echo "Usage: $0 <path_to_forward_file> <path_to_reverse_file> num_split"
	exit 1
elif ! test -f $1
then
	echo $1 is not a valid path
elif ! test -f $2
then
	echo $2 is not a valid path
elif (! grep '^[0-9]*$' <<< $3 > /dev/null) || test $3 -le 0
then
	echo The number of split should be a numeric value greater than 0
fi

# parse command line arguments
forward_path=$1
reverse_path=$2
num_split=$3

#check line number
total_F=`wc -l < "$forward_path"`
total_R=`wc -l < "$reverse_path"`
if test ${total_F} -ne ${total_R}
then
	echo The forward read and the reverse read have different number of lines
	exit 1
fi

# generate file name
forward_name=`echo ${forward_path} | cut -d "." -f1` 
reverse_name=`echo ${reverse_path} | cut -d "." -f1 `

# prepare the directories
mkdir -p Data Intermediate Output/Logfile
mv "$forward_path" "$reverse_path" Data

# calculate the number of lines in each file
num_record=`expr $total_F / 4`
if (( $num_split == 1)) || (( $num_record % 4 == 0))
then
	(( line = ($num_record / $num_split) * 4))
else
	((line = 4 * (($num_record / $num_split) + 1)))
fi

# split the files
split -l${line} -d Data/${forward_name}* "${forward_name}"
split -l${line} -d Data/${reverse_name}* "${reverse_name}"
mv ${forward_name}* Intermediate
mv ${reverse_name}* Intermediate

# make the submission file
cat > submission <<\EOF
#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=6:00:00
#PBS -o output.out
#BPS -e output.err

module load intel/15.0.2
module load python/3.5.1

cd $PBS_O_WORKDIR
EOF
((max=$num_split - 1))
for i in `seq -w 00 $max`
do
	echo "(sh Script/Parser/parse.sh Intermediate/${forward_name}$i Intermediate/${reverse_name}$i Output/LogFile/log_file${i}.txt >> Intermediate/parsed_sequences.txt")& >> submission
done
	echo 'wait&' >> submission

echo Finished

