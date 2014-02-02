#!/usr/bin/perl
#
# This script create an empty unsio sqlite3 database
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
my $db; #     = "/pil/programs/DB/simulation.dbl";  # Marseille database
my $help   = 0; 

GetOptions( 'db=s'        => \$db,
	    'help!'       => \$help,
) or die "Incorrect usage!\n";

# check for help
if( $help  ) {
    usage();
    exit();
}
# get input parameters not named

if ( @ARGV == 1 ) {
    chomp ($db  =$ARGV[0]);
    shift;
} else {
# check mandatory parameters
    die "Option --db not specified.\n" unless defined($db);
}

Sqlite->createUnsioDb($db);
printf STDERR "\nUnsio sqlite3 database [$db] has been created ! \n\n";

# ------------------------------------------------------------
sub usage {
  
  printf STDERR <<EOF;

  ========================================================
  Create an empty unsio sqlite3 database
  ========================================================

 Usage   : $exe databasename
 Example : $exe --db=mydb.sql3

EOF
Tools->printUnsioInfo();
}

