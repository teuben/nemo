#! /bin/sh
# This test is so simple there is no associated .mk file

# The following if/fi block gives you flexibility in choosing the pathname
# for make(1).
if
	test "$1" != ""		# the make pathaname was passed as an argument
then
	MAKE_THIS_SCRIPT="$1" ; export MAKE_THIS_SCRIPT	# so use the pathname
	if
		test ! -x "$MAKE_THIS_SCRIPT" # if the pathname is wrong, stop right now
	then
		echo "$MAKE_THIS_SCRIPT does not exist as an executable file."
		exit 2
	fi
elif
	test ! -z "$MAKE"		# otherwise, check for environment variable
then
	if
		test ! -x "$MAKE"	# if the pathname is wrong, stop right now
	then
		echo "You have MAKE defined as an environment variable, but"
		echo "it does not exist as an executable file."
		exit 2
	else
		echo "Warning: since you have set MAKE explicitly in your environment,"
		echo "This script will say it is defined as a macro."
		echo "To test whether make(1) defines the macro on its own, take MAKE out of your"
		echo "environment and pass the name on the command line to this script, instead."
	MAKE_THIS_SCRIPT=${MAKE} ; export MAKE_THIS_SCRIPT
	fi
else
	# no variable, so let the user's path decide
	MAKE_THIS_SCRIPT=/bin/make ; export MAKE_THIS_SCRIPT
fi

if
	$MAKE_THIS_SCRIPT -pnf - < /dev/null 2>&1 | grep MAKEFLAGS > /dev/null 2>&1
then
	echo 'This version of make defines the MAKEFLAGS macro'
else
	echo 'This version of make does not define the MAKEFLAGS macro'
fi

if
	$MAKE_THIS_SCRIPT -pnf - < /dev/null 2>&1 | grep MFLAGS > /dev/null 2>&1
then
	echo 'This version of make defines the MFLAGS macro'
else
	echo 'This version of make does not define the MFLAGS macro'
fi

if
	$MAKE_THIS_SCRIPT -pnf - < /dev/null 2>&1 | grep 'MAKE *=' > /dev/null 2>&1
then
	echo 'This version of make defines the MAKE macro'
	exit 0
else
	echo 'This version of make does not define the MAKE macro'
	exit 1
fi
