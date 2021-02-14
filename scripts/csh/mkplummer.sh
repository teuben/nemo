#!/bin/bash

trap 'echo "ABORTING";exit' 1 2 15

VERSION=0.1

nbody=100
out=plummer.dat

while getopts n:o:vh i
do
	case $i in
		 h) 	echo "mkplummer demo v$VERSION for NEMO"
				echo "Usage:"
				echo "mkplummer.sh [-n nbody] [-o file] [-h] [-v]"
				echo "n : number of bodies"
				echo "f : filename"
				echo "h : print this help"
				echo "v : print version number"
				exit ;;
		 v) echo "mkplummer.sh  Version $VERSION for NEMO" ; exit ;;
		 n) nbody=$OPTARG ;;
		 o) out=$OPTARG ;;
	esac
done

if [[ -e $out ]]; then
  echo Removing existing $out
  rm -f $out
fi

echo Writing $out with $nbody bodies
mkplummer out=$out nbody=$nbody
