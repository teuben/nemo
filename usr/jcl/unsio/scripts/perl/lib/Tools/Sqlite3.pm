# -*- Mode: Perl -*-
#
# Perl package to manage sqlite3 operations 
#
# ============================================================================
# Copyright Jean-Charles LAMBERT - -2013                                   
# e-mail:   Jean-Charles.Lambert@oamp.fr                                      
# address:  Dynamique des galaxies                                            
#           Laboratoire d'Astrophysique de Marseille                          
#           Pole de l'Etoile, site de Chateau-Gombert                         
#           38, rue Frederic Joliot-Curie                                     
#           13388 Marseille cedex 13 France                                   
#           CNRS U.M.R 6110                                                   
# ============================================================================
#
package Sqlite3::Sqlite;

use strict;
use IO::File;
use Tools::Tools;
use DBI;

my $db="/pil/programs/DB/simulation.dbl";  # Marseille database
my $paramfile="$ENV{HOME}/.unsio";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create unsio sqlite3 empty table
sub Sqlite::createUnsioDb()
{
    my $class          = shift;
    my $db             = shift;

# open database
    my %attr =( PrintWarn=>0,PrintError=>0,RaiseError=>0,Taint=>0);

    my $dbh = DBI->connect( "dbi:SQLite:${db}",,,\%attr ) || die "Cannot connect: $DBI::errstr";
    if ( $dbh->err()) {
      printf STDERR "\nUnable to open DATABASE ==> [$db]\n\n";
      exit();
    }

    # EPS table
    my $stmt = qq(CREATE TABLE eps (name text unique, gas real, halo real, disk real, bulge real, stars real););
    Sqlite->runSqlCommand($dbh,$stmt);
    # INFO table
    my $stmt = qq(CREATE TABLE info (name text unique, type text, dir text, base text););
    Sqlite->runSqlCommand($dbh,$stmt);
    # NEMORANGE table
    my $stmt = qq(CREATE TABLE nemorange (name text unique, total text, disk text, 
                  bulge text, halo text, halo2 text, gas text, bndry text, stars text););
    Sqlite->runSqlCommand($dbh,$stmt);

$dbh->disconnect();
    
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create unsio sqlite3 empty table
sub Sqlite::runSqlCommand() {
    my $class          = shift;
    my $dbh            = shift; # DB file handler
    my $cmd            = shift;

     my $rv = $dbh->do($cmd);
    if($rv < 0){
	print $DBI::errstr;
    } else {
	print "Sql command successfully commited\n";
    }
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Return valid unsio sqlite3 database
sub Sqlite::getValidDb
{
    my $class          = shift;
    my $dbparam        = shift;

    if (defined($dbparam)) {
	$db = $dbparam;
    } else {
	# parse parameter file
	my (%param,$status)=File->parseParameters($paramfile);
	#Tools->printKeyValueHash(\%param); # display parameters
	
	if (defined ($param{dbname}) ) {
	    $db = $param{dbname};
	}
    }
    return $db;
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Return Sim information
#
sub Sqlite::getAllRef {
  # parameters
    my $class          = shift;
    my $simname        = shift;
    my $dbparam        = shift;

    my $db = Sqlite->getValidDb($dbparam);
    die "Unable to open sqlite3 database [$db] ...." if ( ! -f $db );

# open database
    my %attr =( PrintWarn=>0,PrintError=>0,RaiseError=>0,Taint=>0);

    my $dbh = DBI->connect( "dbi:SQLite:${db}",,,\%attr ) || die "Cannot connect: $DBI::errstr";
    if ( $dbh->err()) {
      printf STDERR "\nUnable to open DATABASE ==> [$db]\n\n";
      exit();
    }
    # try Select request
    #my $sql="select * from info where name=='$simname'";
    my $sql="select * from info where name like '$simname' ";
    my $status=$dbh->selectall_hashref($sql,'name');

    if (0) {
      printf "-------------\n";
      foreach my $sim ( $status ) {
	print "sim dir ==> $sim->{dir}\n";
      }
      printf "-------------\n";
    }

    #if (  defined(%$status ) ) {
    if (  %$status  ) {
      return %$status;
    } else {
      return;
    }
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Return Sim information
#
sub Sqlite::getSimInfo {
  # parameters
    my $class          = shift;
    my $simname        = shift;
    my $dbparam        = shift;

    my $db = Sqlite->getValidDb($dbparam);

    die "Unable to open sqlite3 database [$db] ...." if ( ! -f $db );

# open database
    my %attr =( PrintWarn=>0,PrintError=>0,RaiseError=>0,Taint=>0);

    my $dbh = DBI->connect( "dbi:SQLite:${db}",,,\%attr ) || die "Cannot connect: $DBI::errstr";
    if ( $dbh->err()) {
      printf STDERR "\nUnable to open DATABASE ==> [$db]\n\n";
      exit();
    }
    # try Select request
    my $sql="select * from info where name=='$simname'";
    #my $sql="select * from info";
    my $status=$dbh->selectrow_hashref($sql);

    if (0) {
      printf "-------------\n";
      foreach my $sim ( $status ) {
	print "sim dir ==> $sim->{dir}\n";
      }
      printf "-------------\n";
    }

    #if (  defined(%$status ) ) {
    if (  %$status  ) {
      return %$status;
    } else {
      return;
    }
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Return Sim Eps
#
sub Sqlite::getSimEps {
  # parameters
    my $class          = shift;
    my $simname        = shift;
    my $dbparam        = shift;

    my $db = Sqlite->getValidDb($dbparam);
    die "Unable to open sqlite3 database [$db] ...." if ( ! -f $db );

# open database
    my %attr =( PrintWarn=>0,PrintError=>0,RaiseError=>0,Taint=>0);

    my $dbh = DBI->connect( "dbi:SQLite:${db}",,,\%attr ) || die "Cannot connect: $DBI::errstr";
    if ( $dbh->err()) {
      printf STDERR "\nUnable to open DATABASE ==> [$db]\n\n";
      exit();
    }
    # try Select request
    my $sql="select * from eps where name=='$simname'";
    my $status=$dbh->selectrow_hashref($sql);

    if (0) {
      printf "-------------\n";
      foreach my $sim ( $status ) {
	print "sim dir ==> $sim->{dir}\n";
      }
      printf "-------------\n";
    }

    #if (  defined(%$status ) ) {
    if (  $status ) {
      return %$status;
    } else {
      return;
    }
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Return Sim Nemo Range
#
sub Sqlite::getSimNemoRange {
  # parameters
    my $class          = shift;
    my $simname        = shift;
    my $dbparam        = shift;

    my $db = Sqlite->getValidDb($dbparam);
    die "Unable to open sqlite3 database [$db] ...." if ( ! -f $db );

# open database
    my %attr =( PrintWarn=>0,PrintError=>0,RaiseError=>0,Taint=>0);

    my $dbh = DBI->connect( "dbi:SQLite:${db}",,,\%attr ) || die "Cannot connect: $DBI::errstr";
    if ( $dbh->err()) {
      printf STDERR "\nUnable to open DATABASE ==> [$db]\n\n";
      exit();
    }
    # try Select request
    my $sql="select * from nemorange where name=='$simname'";
    my $status=$dbh->selectrow_hashref($sql);

    if (0) {
      printf "-------------\n";
      foreach my $sim ( $status ) {
	print "sim dir ==> $sim->{dir}\n";
      }
      printf "-------------\n";
    }

    #if (  defined(%$status ) ) {
    if (  $status ) {
      return %$status;
    } else {
      return;
    }
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Return Snasphots list
#
sub Sqlite::getSnapshotsList {
  # parameters
    my $class          = shift;
    my $simname        = shift;
    my $dbparam        = shift;
    my $nolog          = shift;

    my $db = Sqlite->getValidDb($dbparam);

    my %res=Sqlite->getSimInfo($simname);

    #if (defined(%res)) {
    if (%res) {
      my $simtype=$res{type};
      #printf STDERR "NOLOG[$nolog]\n";
      # call the right method
      $simtype->getSnapshotsList(\%res,$nolog);
    }
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# check if simulation exist and abort if requested                           
# return simulation's info structure
sub Sqlite::checkSimAbort {
  # parameters
  my $class          = shift;
  my $name           = shift;
  my $abort          = shift;
  
  # collect information from database
  my %res=Sqlite->getSimInfo($name);
  my $status=1;
  #if ( !defined(%res)) {
  if ( !(%res)) {
      $status=0;
    printf STDERR <<EOF;

  Sorry but the simulation $name does not exist in
  [$db] database.....

EOF
    if ($abort) {
      exit();
    }
  }
   return ($status,%res);
}
