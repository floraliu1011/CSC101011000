set -e
if test $# -lt 2
then
	echo "Usage: $0 <path_to_forward_read> <path_to_reverse_read> [optional: <path_to_log_file>"
	exit 1
fi
forward_path=$1
reverse_path=$2
if ! test -z t3
then
	logging_path=$3
fi

#The default script implemented in python
python parser.py ${forward_path} ${reverse_path} ${logging_path}

# if you want to change the program invoked, add you own implementation below and comment the line above
# your script should take in the path to forward read and the path to reverse read 
# and print the parsed sequences into stdout


