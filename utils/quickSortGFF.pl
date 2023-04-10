#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
#use lib '/home/sguizard/work/scripts/seqTools/Module_GFF';
use GFFtools::GFF;
use Bio::SeqIO;
use Getopt::Long;
use Term::ANSIColor;

#> Setting Parameters
##> Define outputs colors
print STDOUT color 'blue';
print STDERR color 'red';

##> Define options
my %config;
GetOptions (\%config,
            'gff=s',
            'sort=i',
            'output_file=s',
            'verbose',
            'help');

##> Print USAGE if --help
if ($config{help}) {printUsage(1);}

##> Check if no mandatory parameter is missing and set gff basename if needed
if (!exists $config{gff})         {printError ("gff option is MANDATORY ! \n", 0);              printUsage();}
if (!exists $config{sort})        {printError ("sort  option is MANDATORY ! \n", 0); printUsage();}
if (!exists $config{output_file}) {printError ("output_file option is MANDATORY ! \n", 0);      printUsage();}
##> Setting Global Variables
$| = 1;
#> Setting parameters


print "Loading gff ... " if $config{verbose};
my $gff=GFF->new(file => $config{gff});
print "OK\n" if $config{verbose};

while (my $seq = $gff->next_seq) {
    if ($config{sort} == 1) {$seq->sort_annot("start","asc");}
    if ($config{sort} == 2) {$seq->sort_annot("start","desc");}
    if ($config{sort} == 3) {$seq->sort_annot("end","asc");}
    if ($config{sort} == 4) {$seq->sort_annot("end","desc");}
}

$gff->write2file($config{output_file});

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
"USAGE : quickSortGFF.pl -g /path/to/gff -s 1 -o sorted_gff.gff3
    
Options : 
    -g | gff                gff file name
    -s | sort               sort by :
                                1 -> start ascending
                                2 -> start descending
                                3 -> end ascending
                                4 -> end descending
    -o | output_file        genome size in bp
    -v | verbose            MORE text dude !!!!
    -h | help               You already know it ...\n\n
    
Sort annotations of a GFF file. 
";
    exit;
}