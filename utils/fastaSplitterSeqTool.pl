#!/usr/bin/perl -w
use strict;
use diagnostics;
use warnings;
use Getopt::Long;
use Term::ANSIColor;
use Data::Dumper;


#> Setting Parameters

##> Define outputs colors
print STDOUT color 'blue';
print STDERR color 'red';


##> Define options
my %config;
GetOptions (\%config,
            'fasta=s',
            'help|h',
            'verbose|v');


##> Print USAGE if --help
if ($config{help}) {printUsage(1);}


##> Check if gff file exists, if no mandatory parameter are missing
if (!exists $config{fasta}){printError ("fasta options are MANDATORY ! \n", 0); printUsage(1);}

##> Main
print "Open fasta file\n" if $config{verbose}; 
open(FASTA, "<$config{fasta}") or printError("Cannot open $config{fasta}", 1);

print "Splitting file\n" if $config{verbose};
while (<FASTA>) {
    chomp;
    if (/>/) {
        no warnings;
        close OUT if (tell(OUT) != -1);
        open(OUT, ">$_.fa") or printError("Cannot open $_.fa", 1);
    }
    print OUT $_."\n";
}

close FASTA;
print "Finish ! \n" if $config{verbose};

###########################################################################
################################ Fonctions ################################
###########################################################################
sub printError{
    my $string = shift;
    my $exit = shift;
    
    print STDERR $string;
    exit if $exit;
}

###########################################################################
sub printUsage{
    my $exit = shift;
    
    print STDOUT
"USAGE : fastaSplitterSeqTool.pl -f a_huge_fasta_file.fa
    
Options :
    -fasta   | -f            fasta file to split (MANDATORY)
    -verbose | -v            MORE text dude !!!!
    -help    | -h            You already know it ...
    
Split a multiple fasta file and ans keep sequence header case for files title.


";

    exit if $exit;
}