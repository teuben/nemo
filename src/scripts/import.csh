#! /bin/csh  -f
#
#	import: 
#           Imports a tar file into an existing ABSOLUTE directory tree
#           Does a number of sanity checks, but is NOT fool proof yet.
#
#	21-mar-90	Created			Peter Teuben
#	29-mar-90	verbosity #files remaining	PJT
#	 7-jun-94       allow 'mkdir -p'		pjt
#	 7-feb-95       also allow tar to be .Z or .gz  pjt
#       13-feb-97       allow directories, and changed
#                       searching directories/files
#       25-mar-97       fixed typo
#	12-apr-97	fixed type in chmod, added optional timestamp update
#	21-nov-98	replaced 'q' with 's' (skip) option	PJT
#

echo "IMPORT: Version 21-nov-98 PJT"

onintr int1_done
if ($#argv != 4) then
    echo "Usage: $0 source destination program mode"
    echo "Copies/Moves all files from a tar-file to a remote dir"
    echo "Directory structure on remote-dir must already exist"
    echo "  <source>:        tar file (can also be .Z or .gz or a directory)"
    echo "  <destination>:   (root) directory (ABSOLUTE)"
    echo "  <program>:       program to use for copy"
    echo "  <mode>:          mode to copy/move the stuff"
    echo "      mode=0       all automatic (untested)"
    echo "      mode=1       plus user menu"
    exit 0
endif
set tar=$1
set remdir=$2
set program="$3"
set mode=$4

set mkdir="mkdir -p"

# echo "IMPORT> tar=$tar dir=$remdir program=$program mode=$mode"

if (! -e $tar) then
  echo "IMPORT Tar file $tar does not exist?"
  exit 0
endif

if (! -d $remdir) then
    echo "Directory root $remdir to import files into does not exist?"
    exit 0
endif

if (-d $tar) then
  set tmp=.
  cd $tar
else
  #       explode files in a temp directory
  set tmp=tmp$$
  $mkdir $tmp
  cd $tmp
  if ($tar:e == gz || $tar:e == Z || $tar:e == tgz) then
    echo "Uncompressing and exploding the tarfile, hang on..."
    zcat ../$tar | tar xf -
  else 
    echo "Exploding the tarfile, hang on..."
    tar xf ../$tar
  endif
  onintr int2_done
endif


#       get list of files to work on, but do this for each directory
#       to prevent the shell from exhausting space
set dirs=(`find . -type d -print`)
set ndirs=$#dirs
foreach dir ($dirs)
  @ ndirs--
  set files=(`echo $dir/*`)
  set nfiles=$#files
  if ($status) then
    echo Skipping seemingly empty directory $dir
    continue
  endif
foreach file ($files)
  @ nfiles--
  if (-d $file) then
    echo skipping directory $file
    continue
  endif
  echo -n "[$ndirs/$nfiles] $file :"
  set rdir=$remdir/$file:h
  if (! -d $rdir) then
    echo ""
    echo "Warning, directory $rdir does not exist. Attempt to create it"
    $mkdir $rdir
    if ($status) continue
    echo "    Created...cross your fingers"
  endif
  if ($mode) then
     cmp -s $file $rdir/$file:t 
     if ($status == 0) then
        echo -n " Identical files - skipped"
        echo ""
        goto loop1_done
     else
        if (! -e $rdir/$file:t) then
            echo " Warning: Destination does not exist"
        else
            echo ""
        endif
     endif
loop1:
     echo -n "Give option: [c]mp [d]iff [l]s [t]ype [q]uit [r]m [s]kip [x]ecute [h]elp:"
     set ans=$<
     switch ($ans)
        case c*:
            cmp $file $rdir/$file:t     
            breaksw
        case d*:
            diff $file $rdir/$file:t | more
            breaksw
        case h*:
        case l*:
            ls -l $file $rdir/$file:t     
            breaksw
	case q*:
            goto int1_done
	case r*:
            rm $file
            goto loop1_done
        case s*:
            goto loop1_done
        case t*:
            more $file 
            breaksw
        case x*:
	    # echo "IMPORT>> $program $file $dir"
            $program $file $rdir
            chmod g+w $rdir/$file:t
            goto loop1_done
        default:
            echo "Illegal option, try again:"
     endsw
     goto loop1
loop1_done:
  else
     cmp -s $file $rdir/$file:t 
     if ($status == 0) then
        echo "Identical files .. skipped"
     else
        $program $file $rdir
        chmod g+w $rdir/$file:t
     endif
  endif
end
end
int2_done:

cd ..
int1_done:

if ("$tmp" != ".") then
  echo "Now deleting temporary directory $tmp..."
  rm -rf $tmp
  echo "and renaming the imported tar file to $tar-"
  mv $tar $tar-
else
  echo ""
  echo All done with directory $tar
endif

echo -n "Update timestamp for adm/chkfile in $NEMO [y/n]? "
set ans=$<
if ($ans == y) then
    cd $NEMO
    make new_time
endif

