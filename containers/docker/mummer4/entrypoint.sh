#!/bin/bash

function usage() {
    # simple entrypoint usage
	echo "Usage: $0 [-h] [-run | -parse] <args ...>"
	echo " "
	echo "Entrypoint for mummer4. use the following commands to select a script."
	echo ""
	echo "Options: "
	echo "          -h                                   Show this help message "
	echo "          -run                                 execute runner script"
	echo "          -shell                               initialize interactive bash shell"
	echo "          <args ...>                           specified arguments for each mode. (Can pass --help)"
	echo " "
}

RUN_MP="simple_ava.sh"

if [[ $# -eq 0 ]]; then
	usage
	echo "Error: Please provide an argument. To view usage for a specific mode, call only that mode."
	exit 1
fi

case $1 in
	"-h")
		usage
		exit 0
		;;
	"-run")
		if [[ $2 == "-shell" ]]; then
			echo "ERROR: Can only specify one mode. Check your command and try again!"
			exit 1
		fi
		shift
		exec /bin/bash "$RUN_MP" "$@"
		;;
    "-shell")
		exec /bin/bash
		;;
	*)
		echo "ERROR: invalid argument! Please only use either '-run' or '-parse'"
		exit 1
		;;
esac

exit 0