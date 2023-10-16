#! /usr/bin/env bash
#
#     can also try:    /bin/bash
#

usage() {
    echo "Usage: $0 [-a AAA] [-b BBB]"
    exit 1;
}

# if no arguments are allowed, give usage and exit
[ $# -eq 0 ] && usage

while getopts ":ha:b:" arg; do
  case $arg in
    a) # Specify aaa value.
	echo "aaa is ${OPTARG}"
	aaa=${OPTARG}	
      ;;
    b) # Specify bbb value
	bbb=${OPTARG}
      ;;
    h | *) # Display help.
      usage
      exit 0
      ;;
  esac
done

echo "Computing with aaa=${aaa} and bbb=${bbb}"
