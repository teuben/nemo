#! /usr/bin/bash
#
#  If you have "about equal CPU" tasks that are normally single CPU programs,
#  running them with this simply shell construct can be useful to speed up
#  a series
#


echo Sleeping 2
sleep 2 &

echo Sleeping 4
sleep 4 &

echo Sleeping 6
sleep 6 &

echo Waiting until all are done....
wait

echo All done
