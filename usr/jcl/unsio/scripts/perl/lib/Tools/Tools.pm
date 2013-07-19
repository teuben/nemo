# -*- Mode: Perl -*-
#
# Perl package with useful functions
#
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
package Tools::Tools;

use strict;
use IO::File;

# - - - - - - - - - - - - - - - - - - - - - - - - - -
#
sub Tools::printUnsioInfo
{
    my $class  = shift;

    printf STDERR <<EOF;

  To get information about UNSIO, please visit : 
  -> http://projets.lam.fr/projects/unsio

  To get information about Unsio Sqlite3 database, please visit :
  -> http://projets.lam.fr/projects/unsio/wiki/Sqlite3Db

EOF
}
# - - - - - - - - - - - - - - - - - - - - - - - - - -
#
sub File::readFile
{
    my $class  = shift;
    my $HANDLE = shift;
    my $line;

    #printf STDERR "In Tools2 [$HANDLE]\n";
    while ( <$HANDLE> ) {
	chomp;
	#printf STDERR ">>> %s\n", $_;
	return (1,$_); 
    }
    return 0;
}

# - - - - - - - - - - - - - - - - - - - - - - - - - -
# parse a file with input parameters based on pair "key = value"
# and fill up the result into an hash table
sub File::parseParameters
{
    my $class  = shift;
    my $file    = shift;  # input file

    my %param;  # hash structure to store parameter file

    my $FD =  IO::File->new();

    if 	(! $FD->open("<$file")) {
	warn "Unable to open for reading [$file]";
	return (%param,0);
    }

    my $ret=1;

    while ($ret) {
	my $line;
	($ret,$line)  = File->readFile($FD);
	if ($ret) {
	    my $line2 = (split(/#/,$line))[0]; # skip commented line
	    if (defined($line2)) {
		my ($key, $value) = split(/=/,$line2);
		if ((defined $key) && (defined $value)) {
		    $key   =~ s/ +//g;   # remove all blank characters
		    $value =~ s/^\s+//;  # remove leading  blank characters
		    $value =~ s/\s+$//;  # remove trailing blank characters
		    #printf STDERR "[$key]/[$value]\n";
		    $param{$key} = $value;
		}
	    }
	    #my $simname="$tab[0]";  # simulation name
	}
    }
    return (%param,1);
}
# - - - - - - - - - - - - - - - - - - - - - - - - - -
#
sub Tools::printKeyValueHash
{
    my $class = shift;
    #my %hash  = %{(shift)}; #% hash table given in parameter (copy in local variabl)
    my $hash  = shift; #% hash table given in parameter (pass by reference)

    #foreach my $kk ( keys((%hash)) ) { #
    foreach my $kk ( keys(($hash)) ) { #
	printf STDERR "Key/Value = [$kk/$hash->{$kk}]\n";
    }
}

# - - - - - - - - - - - - - - - - - - - - - - - - - -
#
sub Tools::run
{
    my $class      = shift;
    my $cmd        = shift;
    my $notverbose = shift;
    
    if (! defined($notverbose)) {
	print STDERR "Running:\n[$cmd]\n";
    }
    system($cmd);
}
# - - - - - - - - - - - - - - - - - - - - - - - - - -
# Guess CPU number 
# and create hostmpich
sub guess_cpu_number
{

    #system("echo `hostname|cut -d. -f1`:$ncpu > hostmpich");
    my $host=Tools->getMyShortHostname();
    my $nproc=Tools->getCoresNumber();

    my $hmp="hostmpich";
    my $HMP = IO::File->new("> $hmp ") ||
        die "Unable to open for writing [$hmp]";
    for (my $i=0; $i<${nproc}; ${i}++) {
	printf $HMP "${host}\n";
    }
    $HMP->close();
        
    printf "#CPU detected = $nproc\n";
    return $nproc;
}
# - - - - - - - - - - - - - - - - - - - - - - - - - -
# 
sub Tools::getCoresNumber {
  # parameters
  my $class = shift;

  my $cpu=0;
  my $FD;
  if ( open($FD,"cat /proc/cpuinfo | grep processor|wc -l|")) {
    chomp(my $x=<$FD>);
    if ( defined $x ) {
	$cpu=$x;
    }
    close ($FD);
  } else {
    printf STDERR "In Tools::getCoresNumber, unable to detect #cores..\n";
  }
  return $cpu;
}
# - - - - - - - - - - - - - - - - - - - - - - - - - -
1;  # Allways return 1

#
