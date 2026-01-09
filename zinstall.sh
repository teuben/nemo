#!/bin/sh
#
#  graphical aid, using zenity, to configure and build NEMO
#  https://gitlab.gnome.org/GNOME/zenity
#  see also https://yad-guide.ingk.se/


if ! [ -x "$(command -v zenity)" ]; then
  echo 'Error: zenity is not installed.' >&2
  exit 1
fi



options=$(zenity --list --checklist \
       --title="Choose the configure and build options for NEMO" \
       --separator=, \
       --column="Option"  --column "what"  --column="Description" \
       FALSE   enable-debug      "Debugging enables" \
       FALSE   disable-fortran   "Disable footran" \
       FALSE   with-csh          "no csh" \
       FALSE   disable-shared    " on a mac?" \
       FALSE   enable-clang      "clang compiler instead of gcc" \
       FALSE   enable-single     "single precision" \
       TRUE    build1    "do this" \
       TRUE    build2    "do this" \
       TRUE    build3    "do this" \
       FALSE   build4    "do this" \
       TRUE    check     "do this" \
       TRUE    bench5    "do this")


echo $options

run=my_zinstall.sh


echo "#  $run : an install file for NEMO as created by $0"  > $run
echo "./configure"  >> $run
echo "make build1"  >> $run
echo "make build2"  >> $run
echo "make build3"  >> $run
echo "make build4"  >> $run
echo "make check"   >> $run
echo "make bench5"  >> $run

zenity --text-info \
       --title="Order of commands for installing NEMO" \
       --filename=$run \
       --checkbox="Go ahead and run this install"

case $? in
    0) echo "Starting the installation"
       ;;
    1) echo "Stop"
       exit 0
       ;;
   -1)
       echo "Unexpected error"
       ;;
esac

chmod +x $run
echo "Not quite ready to do this automagically...."
echo "./$run"
