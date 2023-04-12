#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use Getopt::Long;

###> Loading Options
my %config;
GetOptions (\%config,
            'fa=s',
            'out',
            'help');

if ($config{help})        {printUsage(1);}
if (!exists $config{fa})  {printError("fa option is mandatory ! \n", 1);}

open(my $fh_in, "<", $config{fa}); 

my $tfa; 
my $fh_out; 

if ($config{out}) {
    $tfa = $config{fa}; 
    $tfa =~ s/(fa|fasta)$/tfa/; 
    open($fh_out, ">", $tfa);
}

while (<$fh_in>) {
    if ($. == 1) {
        s/\n$/\t/; 
        printOut($_);
        next;
    } 
    if (/^>/) {
        $_ = "\n$_"; 
        s/\n$/\t/g; 
        printOut($_);
        next;
    } 
    chomp;
    printOut($_);
}

close $fh_in;
close $fh_out if ($config{out});


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
sub printOut{
    my $string = shift;
    
    if ($config{out}) {
        print $fh_out $string;
    } else {
        print STDOUT $string;
    }
}

###########################################################################
sub printUsage{
    my $exit = shift;
    
    print STDOUT
"USAGE : fa2tfa.pl -fa seq.fa
    
Options : 
    -fa     A fasta file (string) (Mandatory)
    -out    Save data in .tfa file
    -help   You already know it ...

Convert regular fasta in tabular fasta (col 1: seq id; col 2: seq). 

";
    exit if $exit;
}