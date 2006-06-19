#!/usr/bin/perl
# ============================================================================|
# Copyright Jean-Charles LAMBERT - 2004-2006                                  |
# e-mail:   Jean-Charles.Lambert@oamp.fr                                      |
# address:  Dynamique des galaxies                                            |
#           Laboratoire d'Astrophysique de Marseille                          |
#           2, place Le Verrier                                               |
#           13248 Marseille Cedex 4, France                                   |
#           CNRS U.M.R 6110                                                   |
# ============================================================================|
# See the complete license in LICENSE and/or "http:#www.cecill.info".         |
# ============================================================================|
#                                                                             |
#                   Create DivX movie from a set of PNG files                 |
#                                                                             |
# ============================================================================|

use strict;

# get exe basename
my $exe=(split(/\//,$0))[-1];

# Check parameters
if (@ARGV < 5 ) {
  printUsage();
  exit();
}

# get parameters
chomp (my $dir   = $ARGV[0]);
shift (@ARGV);
chomp (my $base  = $ARGV[0]);
shift (@ARGV);
chomp (my $fps   = $ARGV[0]);
shift (@ARGV);
chomp (my $movie = $ARGV[0]);
shift (@ARGV);
chomp (my $codec = $ARGV[0]);
shift (@ARGV);
chomp (my $bitfactor = $ARGV[0]);
shift (@ARGV);
if ($bitfactor == "") {
  $bitfactor = 1;
}


# get parameters
my %info=getInfoCodec();
if ( ! $info{$codec} ) {
  printf STDERR "\nCodec [$codec] does not exist, exiting....\n\n";
  exit();
}

createMovie(\%info, $dir, $base, $fps, $movie, $codec, $bitfactor);

# -----------------------------------------------------------
sub createMovie() {
  my $info = shift;
  my $dir  = shift;
  my $base = shift;
  my $fps  = shift;
  my $movie= shift;
  my $codec= shift;
  my $bitfactor = shift;
 
  my $optbitrate;
  my $cpt=0;
  foreach ( <${dir}/${base}*.png> ) {
    $cpt++;
    if ($cpt==1) {
      my ($W,$H) = getWH($_);
      printf STDERR "[$_] $W  X  $H\n";
      $optbitrate = int($bitfactor * 50 * 25 * $W * $H / 256);
    }
  }
  # mencoder (pass1 and pass2)
  my $menc1="mencoder -ovc lavc -lavcopts vcodec=${codec}:vpass=1:vbitrate=${optbitrate}:$info->{$codec}->{options}";
  my $menc2="mencoder -ovc lavc -lavcopts vcodec=${codec}:vpass=2:vbitrate=${optbitrate}:$info->{$codec}->{options}";
  # mencoder frame options
  my $mencfps="-mf type=png:fps=${fps} -nosound";
  # mencoder output (pass1 and pass2)
  my $menco1="-o /dev/null \"mf://${dir}/${base}*.png\"";
  my $menco2="-o ${movie} \"mf://${dir}/${base}*.png\"";

  if (-f "divx2pass.log" ) {
    unlink "divx2pass.log";# || warning "Unable to remove [divx2pass.log]\n";
  }

  my $cmd="$menc1 $mencfps $menco1";
  printf STDERR "running [$cmd]\n";
  system($cmd);

  my $cmd="$menc2 $mencfps $menco2";
  printf STDERR "running [$cmd]\n";
  system($cmd);  
}
# -----------------------------------------------------------
sub getWH() {
  my $frame = shift;

  my $FD;
  open($FD,"identify $_ |");
  chomp(my $x=<$FD>);
  $x = (split(/ /,$x))[2];
  my $W=(split(/x/,$x))[0];
  my $H=(split(/x/,$x))[1];


  return $W,$H;
}
# -----------------------------------------------------------
sub getInfoCodec() {

  my %pp = ( 
	    msmpeg4v2 => {
			  options => "mbd=2:keyint=132:vqblur=1.0:cmp=2:subcmp=2:dia=2:mv0:last_pred=3"
			 },
	    mpeg4     => {
			  options => "mbd=2:keyint=132:v4mv:vqmin=3:lumi_mask=0.07:dark_mask=0.2:mpeg_quant:scplx_mask=0.1:tcplx_mask=0.1:naq"
			 }
	   );
  return %pp;
}
# -----------------------------------------------------------
sub printUsage() {

printf STDERR <<EOF;


         ------------- Mencoder front End -------------
           Create DivX movie from a set of PNG files
         ----------------------------------------------


 Usage   : $exe directory   basename_frame fps movie_out_name codec_type [bitfactor]
 Example : $exe /tmp/render frame  25 movie.avi  msmpeg4v2  3.5

bitfactor:
----------
 this value will be multiplied by the optimal bit rate.

available codec_type are:
=========================
 - msmpeg4v2  : Microsoft DivX codec
 - mpeg4      : High quality Dix4/5 codec (might be unreadble on M\$)



Author: Jean-Charles LAMBERT - 2004-2006


EOF
}
