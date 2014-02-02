#!/usr/bin/perl
#
# This script update INFO table from an unsio sqlite3 database
#
#
# ============================================================================
# Copyright Jean-Charles LAMBERT - 2008-2014
#           Centre de donneeS Astrophysiques de Marseille (CeSAM)       
# e-mail:   Jean-Charles.Lambert@lam.fr                                      
# address:  Aix Marseille Universite, CNRS, LAM 
#           Laboratoire d'Astrophysique de Marseille                          
#           Pole de l'Etoile, site de Chateau-Gombert                         
#           38, rue Frederic Joliot-Curie                                     
#           13388 Marseille cedex 13 France                                   
#           CNRS UMR 7326    
# ============================================================================
#
#
#

BEGIN 
{ 
    my $exe=(split(/\//,$0))[-1]; # basename program name
    my $path=$0;                  # absolute path name
    $path =~ s/${exe}//g;         # dirname
    push @INC, "$path/../scripts/perl/lib", "$path/../lib"; # add path to locate modules
}

# load package
use strict;                    # must define everything
use Getopt::Long;
use Tools::Sqlite3;


# -------------------------------------------------------------
# Main program
# get exe basename
my $exe=(split(/\//,$0))[-1];

# setup default input parameter
my $name;
my $db; #     = "/pil/programs/DB/simulation.dbl";  # Marseille database
my $help   = 0; 

GetOptions( 'simname=s'   => \$name,
            'db=s'        => \$db,
	    'help!'       => \$help,
) or die "Incorrect usage!\n";

# check for help
if( $help  ) {
    usage();
    exit();
}
# get input parameters not named
if ( @ARGV >= 1) {
    chomp($name=$ARGV[0]); shift;
    if ( @ARGV == 1 ) {
	chomp ($db  =$ARGV[0]);
	shift;
    }
} else {
# check mandatory parameters
    die "Option --simname not specified.\n" unless defined($name);
}

if ($db eq "") {
    undef $db;
}

my $db = Sqlite->getValidDb($db);

printf STDERR "simname=[$name] db=[$db]\n";

my %res=Sqlite->getSimInfo($name,$db);

if ((%res)) {
  printf "dir = $res{dir}\n";

  if (1) {
    for my $info ( keys %res ) {
      printf " $info = $res{$info}\n";
    }
  }
}

my %res=Sqlite->getSimEps($name,$db);

if ((%res)) {
  if (1) {
    for my $info ( keys %res ) {
      printf " $info = $res{$info}\n";
    }
  }
}

my %res=Sqlite->getSimNemoRange($name,$db);

if ((%res)) {
  if (1) {
    for my $info ( keys %res ) {
      printf " $info = $res{$info}\n";
    }
  }
}

# ------------------------------------------------------------
sub usage {
    my $db = Sqlite->getValidDb($db);
  printf STDERR <<EOF;

  ========================================================
  Return from an unsio sqlite3 database simulation file:
  [$db]
  all the informations belonging to the simulation\'s name
  given in parameter.
  ========================================================

 Usage   : $exe simname [db]
 Example : $exe --simname=sgs019

EOF
Tools->printUnsioInfo();
}

