#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use Getopt::Long;
use Term::ANSIColor;
use Bio::SeqIO;


###> Loading Options
my %config;

GetOptions (\%config, 'fa=s', 'id=s', 'output=s', 'verbose', 'help');

if ($config{help} or (!exists $config{fa} or !exists $config{id} or !exists $config{output})) {
    print color 'blue';
    print
"USAGE : extractByID.pl -fa multi.fasta -id id_list.txt -o output.fa
    
Options : 
    --fa        Multi Fasta file (Mandatory)
    --id        A text file listing id to extract (Mandatory)
    --output    Output to create (Mandatory)
    --verbose   MORE text dude !!!!
    --help      you already know it ...\n\n";
    exit;
}

if (!exists $config{fa})     {printError("fa is mandatory ! \n", 1);}
if (!exists $config{id})     {printError("id is mandatory ! \n", 1);}
if (!exists $config{output}) {printError("output is mandatory ! \n", 1);}

if (!-e $config{fa})         {printError("$config{fa} do not exist ! \n", 1);}
if (!-e $config{id})         {printError("$config{id} do not exist ! \n", 1);}


###> Initialize parameters
$| = 1;
my %seq; 


###> Open multi fasta $config{fa}
my $seqIN  = Bio::SeqIO->new(-format => 'fasta', -file => $config{fa});
while (my $s = $seqIN->next_seq()) {
    $seq{$s->display_id()} = $s; 
}


###> Open output fasta $config{output}
my $seqOUT = Bio::SeqIO->new(-format => 'fasta', -file => ">$config{output}"); 


###> Open id list $config{id}
open(ID, "<$config{id}") or printError("Cannot open $config{id}", 1); 


###> Extracting sequences
while (<ID>) {
    chomp;
    s/>//g;
    
    if (!exists($seq{$_})) {
        printError("Can't find $_", 1);
    }
    
    $seqOUT->write_seq($seq{$_});
}

close ID;


###########################################################################
################################ Fonctions ################################
###########################################################################
sub printError{
    my $string = shift;
    my $exit = shift;
    
    print color 'red';
    print STDERR "ERROR: $string";
    exit if $exit;
}