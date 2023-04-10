#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use Getopt::Long;
use Term::ANSIColor;


#> Setting Parameters
##> Define outputs colors
print STDOUT color 'blue';
print STDERR color 'red';

##> Define options
my %config;
GetOptions (\%config,
            'multi_gff=s',
            'sufixe=s',
            'verbose',
            'help');

##> Print USAGE if --help
if ($config{help})           {printUsage(1);}
if (!exists $config{sufixe}) {$config{sufixe} = "";}


##> Check if no mandatory parameter is missing and set gff basename if needed
if (!exists $config{multi_gff}) {printError ("multi_gff option is MANDATORY ! \n", 0); printUsage();}

##> Setting Global Variables
$| = 1;
#> Setting parameters

open(MULTI, "<", $config{multi_gff}) or printError("Cannot open $config{multi_gff}", 1);

my $ok = 0;

while (<MULTI>) {
    next if /##gff-version 3/;
    if (/##sequence-region (\S+) \d+ \d+/) {
        close OUT if $ok;
        $ok = 1;
        print "Creating file ${1}$config{sufixe}.gff3\n" if $config{verbose};
        open(OUT, ">", "${1}$config{sufixe}.gff3") or printError("Cannot open ${1}$config{sufixe}.gff3", 1);
        print OUT "##gff-version 3\n";
    }
    print OUT if $ok;
}



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

    print STDOUT
"USAGE : quickGffSplitter.pl -m gff.gff3
    
Options : 
    -m | multi_gff              multi gff to process
    -s | sufixe                 suffixe to add to filename
    -v | verbose                MORE text dude !!!!
    -h | help                   You already know it ...
    
split a multi gff into single gff. For more advance GFF splitter use : quickGffSpliterAndFormator.pl

";
    exit;
}

