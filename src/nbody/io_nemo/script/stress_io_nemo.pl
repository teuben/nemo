#!/usr/bin/perl

use strict;

my $NEMOSRC=$ENV{NEMOSRC}; # get $NEMOSRC environment variable
my $OSTYPE=$ENV{OSTYPE};   # get $OSTYPE environment variable
# Check necessary nemo binaries
my @needprg=( "hackcode1", "snapmask", "snapprint", "mkplummer" );
my $FD;
foreach my $i ( @needprg  ) { 
  open($FD,"$NEMOSRC/scripts/need $i |");
  close($FD);
  chomp(my $x=<$FD>);
  printf STDERR "$i [$x]\n";
  if ( $x =~ /NEMOBIN/ ) {
    system("$NEMOSRC/scripts/need -m  $i");
  }
}

# create input snapshot and particle range
my @sim;
my $i=0;
$sim[$i][0] = "/tmp/stress_io_p${i}";
$sim[$i][1] = "500";
$sim[$i][2] = "0:100,200:250";
$i++;
$sim[$i][0] = "/tmp/stress_io_p${i}";
$sim[$i][1] = "800";
$sim[$i][2] = "all";
$i++;
$sim[$i][0] = "/tmp/stress_io_p${i}";
$sim[$i][1] = "1000";
$sim[$i][2] = "0:500,501:999:2";
$i++;
$sim[$i][0] = "/tmp/stress_io_p${i}";
$sim[$i][1] = "5000";
$sim[$i][2] = "501:999,1000:2000,2500:4000:3,4001:4999";
$i++;
$sim[$i][0] = "/tmp/stress_io_p${i}";
$sim[$i][1] = "2000";
$sim[$i][2] = "all";
$i++;

buildModel($i);                   # build models
my $inlist=buildInputList($i);    # build input list

foreach my $p ("stress_io_nemo" , "stress_io_nemo_f") {
  foreach my $f ("double" , "float")  {
    stressIoNemo($inlist,$i,$f,$p);# run stress_io_nemo program
    compareData($i,$f);            # compare result data
    printf STDERR "\n";
  }
}

# --------------------------------------------------------------------
sub compareData {
  my $nmodel = shift;
  my $ftype  = shift;
  
  my $iopath="${NEMOSRC}/nbody/io_nemo/${OSTYPE}/bin";
  my $mask_prog;
  if ( $ftype eq "double" ) {
    $mask_prog="$iopath/snapmask_d";
  }
  else {
    $mask_prog="$iopath/snapmask_s"; 
  }
  printf STDERR "Mask program used: [$mask_prog] with datatype [$ftype] \n";
  #
  for (my $n=0; $n < $nmodel; $n++) {

    # original
    my $in="$sim[$n][0]";
    my $out="$sim[$n][0]_mask";
    `/bin/rm -f $out 2> /dev/null`;
    `($mask_prog $in - select=$sim[$n][2] |\
      snapprint - options=t,m,x,y,z,vx,vy,vz,ax,ay,az,phi > $out) 2>/dev/null`;

    # io_nemo
    my $in2="$sim[$n][0]_out";
    my $out2="$sim[$n][0]_out_mask";
    `/bin/rm -f $out2 2> /dev/null`;
     `(snapprint $in2 options=t,m,x,y,z,vx,vy,vz,ax,ay,az,phi > $out2) 2>/dev/null`;

    # compare
    printf STDERR "Comparing [$out] vs [$out2]...\n";
    my $bad=`diff $out $out2 | wc -l`;
    if ( $bad != 0 || -z $out || -z $out2 ) {
      printf "*ERROR* $out $out2 differ....\n";
      exit();
    }
  }
}
# --------------------------------------------------------------------
sub stressIoNemo {
  my $inlist = shift;
  my $nmodel = shift;
  my $ftype  = shift;
  my $program= shift;
  #
  # get rid off
  for (my $n=0; $n < $nmodel; $n++) {
    my $out="$sim[$n][0]_out";
    `/bin/rm -f $out 2> /dev/null`;
  }  
  my $iopath="${NEMOSRC}/nbody/io_nemo/${OSTYPE}/bin";
  printf STDERR "Running <$program> program [$ftype] datatype.....\n";
  my $cmd="${iopath}/${program} $inlist $ftype 2> /dev/null";
  system($cmd);
}
# --------------------------------------------------------------------
sub buildModel {
  my $nmodel = shift;

  # hackcode1 parameters
  my $tstop=0.3125;
  my $freqout=32;
  my $minor_freqout=32;

  # build models
  for (my $n=0; $n < $nmodel; $n++) {
    my $name="$sim[$n][0]";
    printf STDERR "$sim[$n][0]   $sim[$n][1]   $sim[$n][2]\n";
    #`/bin/rm -f $name 2> /dev/null`;
    if ( ! -f $name) {
      printf STDERR "Building model [$name], please wait ...\n";
      `(mkplummer - ${sim[$n][1]}       |\
        hackcode1 - out=${name} options=mass,phase,phi,acc tstop=$tstop freqout=$freqout minor_freqout=$minor_freqout) 2> /dev/null`
    }
    
  }
}
# --------------------------------------------------------------------
sub buildInputList {
  my $nmodel = shift;
  
  my $FD;
  my $inlist="/tmp/stress_io_input_list";
  printf STDERR "Building input list file....\n";
  open($FD,"> $inlist");
  for (my $n=0; $n < $nmodel; $n++) {
    printf $FD "$sim[$n][0]\n$sim[$n][1]\n$sim[$n][2]\n";
  }
  close($FD);
  return $inlist;
}
# --------------------------------------------------------------------
