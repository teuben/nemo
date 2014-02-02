#!/usr/bin/perl
#
# This script return simulation information from an unsio sqlite3 database
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
use Tools::Tools;
use Tools::Sqlite3;
use DBI;
# -------------------------------------------------------------
# Main program
# get exe basename
my $exe=(split(/\//,$0))[-1];

# setup default input parameter
my $name;
my $type;
my $dir;
my $base;
my $db;
my $help   = 0; 

GetOptions( 'simname=s'   => \$name,
	    'type=s'      => \$type,
            'dir=s'       => \$dir,
            'base=s'      => \$base,
            'db=s'        => \$db,
	    'help!'       => \$help,
) or die "Incorrect usage!\n";

my $x=@ARGV;
#printf STDERR "ARGV = [$x]\n";

# check for help
if( $help || ( @ARGV > 1  && @ARGV < 4)  ) {
    usage();
    exit();
}

# get input parameters not named
if ( @ARGV >= 4) {
    chomp($name=$ARGV[0]); shift;
    chomp($type=$ARGV[0]); shift;
    chomp($dir =$ARGV[0]); shift;
    chomp($base=$ARGV[0]); shift;
    if ( @ARGV == 1 ) {
	chomp ($db  =$ARGV[0]);
	shift;
    }
} else {
# check mandatory parameters
    die "Option --simname not specified.\n" unless defined($name);
    die "Option --type not specified.\n" unless defined($type);
    die "Option --dir not specified.\n" unless defined($dir);
    die "Option --base not specified.\n" unless defined($base);
}

my $db = Sqlite->getValidDb($db);

printf STDERR "simname=[$name] type=[$type] dir=[$dir] base=[$base] db=[$db]\n";

# check if file exist
die "Unsio database [$db] does not exist...." if (! -f $db) ;

# open sqlite3 database
my %attr =( PrintWarn=>0,PrintError=>0,RaiseError=>0,Taint=>0);

my $dbh = DBI->connect( "dbi:SQLite:${db}",,,\%attr ) || die "Cannot connect: $DBI::errstr";

# try to INSERT request
my $sql="insert into info (name,type,dir,base) values ('$name','$type','$dir','$base')";
my $status=$dbh->do($sql);
if ( $dbh->err()) {
  #printf STDERR "Record already exist, trying to update now...\n";
  $sql="update info set name='$name',type='$type',dir='$dir',base='$base' where name='$name'";
  my $status=$dbh->do($sql);
  if ( $dbh->err()) {
    printf STDERR "\nFAILED TO UPDATE\n";
  } else {
    printf STDERR "Record [$name] sucessfully UPDATED.\n";
  }
} else {
   printf STDERR "Record [$name] sucessfully INSERTED\n";
}

# ------------------------------------------------------------
# HOW to use the program
sub usage {

  printf STDERR <<EOF;

      ==========================================================
      Update Info table from an UNSIO SQLite3 database
      ==========================================================

    usage   : $exe   sim_name         sim_type      sim_dir                     sim_basename       [unsio sqlite3 db]
    example : $exe --simname=sgs007 --type=Gadget --dir=/rbdata/sergas/sgs007 --base=SNAPS/sgs007  [--db=simulation.dbl]


  type :   Gadget or Nemo

  Tips: to update from a text file, use the following command:
  cat /pil/programs/DB/sim_info.txt | grep -v # | xargs -l1 ${exe}

  OR

  to update all the sgs??? Gadget simulation located to /radata/sergas

  ls -d /radata/sergas/sgs??? |        \\
       xargs -IXX basename XX |        \\
       xargs -IXX ${exe} \\
       --simname=XX --type=Gadget --dir=/radata/sergas/XX --base=SNAPS/XX --db=/pil/programs/DB/simulation.dbl

EOF
Tools->printUnsioInfo();
}
#
