#	This piece of csh code is supposed to be 'source'd
#	in your shell script, and allows you to use keywords
#	plus defaults.
#
#	Parameters can be supplied in the form of commandline
#	arguments:
#		command KEY1 VALUE1 KEY2 VALUE2 ...
#

##  THIS SECTION IS TO BE SUPPLIED BY THE USER
##  SET KEYS, VALS and HELP shell variables
# set keys=(a b c)
# set vals=(1 2 3)
# set help="Script to create and analyse velocity field"
## END SET

##  Check if anything supplied at all
if ($#argv == 0) then
   echo "Usage: $0 [KEY VALUE ...]"
   echo ""
   @ err=1
else
   @ err=0
endif

##  Loop installing new values for keyword from command line
loop1:
  if ($#argv < 2) goto loop1_done

    @ chk=0
    @ i=0
    foreach key ($keys)
      @ i++
      if ($key == $argv[1]) then
        echo "Resetting key-$i $argv[1] to $argv[2]"
        set vals[$i]=$argv[2]
	@ chk++
      endif
    end
    if ($chk == 0) then
      echo "Unknown key=$argv[1] value=$argv[2]"
      @ err++
    endif

  shift argv
  shift argv
  goto loop1
loop1_done:

## Display and set keywords for script
  echo "Current Values of shell script keywords:"
  @ i=0
  foreach key ($keys)
    @ i++
    echo "$keys[$i] = $vals[$i]"
    set $keys[$i]=$vals[$i]
  end
  if ($#argv > 0) then
    echo ""
    echo "***ERROR: unprocessed remaining arguments: $*"
    exit 0
  endif
  if ($err != 0) then
    echo ""
    echo "$help"
    echo ""
    echo "[You need to supply at least one KEY VALUE pair to start up script]"
    exit 0
  endif
## end set and display keywords

echo "Script $0 running now ...."
#================================================================
#	begin of user customized part of script
#================================================================
