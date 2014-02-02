#!/usr/bin/perl
#
# This script update NEMORANGE table from an unsio sqlite3 database
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
my $total ='';
my $disk  ='';
my $bulge ='';
my $halo  ='';
my $halo2 ='';
my $gas   ='';
my $bndry ='';
my $stars ='';
my $db     = "";  # Marseille database
my $help   = 0; 

GetOptions( 'simname=s'    => \$name,
	    'total:s'      => \$total,
            'disk:s'       => \$disk,
            'bulge:s'      => \$bulge,
            'halo:s'       => \$halo,
            'halo2:s'      => \$halo2,
            'gas:s'        => \$gas,
            'bndry:s'      => \$bndry,
            'stars:s'      => \$stars,
            'db=s'         => \$db,
	    'help!'        => \$help,
) or die "Incorrect usage!\n";

if ( $help ) {
  usage();
  exit();
} else {
    # check mandatory parameters
    die "Option --simname not specified.\n" unless defined($name);
}
my $db = Sqlite->getValidDb($db);

printf STDERR "name=$name,total=$total,disk=$disk,bulge=$bulge,halo=$halo,halo2=$halo2 ,gas=$gas ,bndry=$bndry, stars=$stars [$db]\n";

# check if file exist
die "Unsio database [$db] does not exist...." if (! -f $db) ;

# open database
my %attr =( PrintWarn=>0,PrintError=>0,RaiseError=>0,Taint=>0);

my $dbh = DBI->connect( "dbi:SQLite:${db}",,,\%attr ) || die "Cannot connect: $DBI::errstr";


#exit();

# try to INSERT request
my $sql="insert into nemorange (name,total,disk,bulge,halo,halo2,gas,bndry,stars) values ('$name','$total','$disk','$bulge','$halo','$halo2','$gas','$bndry','$stars')";
my $status=$dbh->do($sql);
if ( $dbh->err()) {
  #printf STDERR "Record already exist, trying to update now...\n";
  $sql="update nemorange set name='$name',total='$total',disk='$disk',bulge='$bulge',halo='$halo',halo2='$halo2',gas='$gas',bndry='$bndry', stars='$stars' where name='$name'";
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
sub usage {

  printf STDERR <<EOF;

      ==========================================================
      Update NEMORANGE table of an UNSIO SQLite3 database
      ==========================================================

       usage   : $exe simname    [total] [disk] [bulge] [halo] [halo2] [gas] [bndry] [stars] [database]
       example : $exe --simname=nde145 --total=0:1932766 --disk=0:99999 --halo=100000:599772 --halo2=599773:1932766

EOF
Tools->printUnsioInfo();
}


#
