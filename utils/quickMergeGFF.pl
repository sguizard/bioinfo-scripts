#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use lib '/home/sguizard/work/scripts/seqTools/Module_GFF_2';
use GFF;
use Annotation;
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
            'verbose',
            'help');

##> Print USAGE if --help
if ($config{help}) {printUsage(1);}

##> Setting Global Variables
$| = 1;
#> Setting parameters

print "Creating genome GFF\n" if $config{verbose};
my $outGff = GFF->new(file => "galgal4_genome.gff3");
my $b      = 1;
my $sumSeq = 0;


foreach my $f (<*.gff3>){
    print "\rProcessing file : $f" if $config{verbose};    
    my $gff= GFF->new(file => $f);
    
    if ($b) {
        $outGff->version(3);
        $outGff->sequence_region->id("galgal4_genome");
        $outGff->sequence_region->start(1);
        
        while ($gff->next_annot) {
            $gff->seq_id("galgal4_genome");
            $outGff->add_annot($gff->get_annot);
        }
        $sumSeq += $gff->sequence_region->end;
        $b--;
    }
    else {
        while ($gff->next_annot) {
            $gff->seq_id("galgal4_genome");
            $gff->start($gff->start + $sumSeq);
            $gff->end($gff->end + $sumSeq);
            $outGff->add_annot($gff->get_annot);
        }
        $sumSeq += $gff->sequence_region->end;
    }
}
print "\nParsing done \n" if $config{verbose};

$outGff->sequence_region->end($sumSeq);

print "Writing file\n" if $config{verbose};
$outGff->write2file;

print "Done\n" if $config{verbose};

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
"USAGE : quickMergeGFF.pl -v
    
Options : 
    -v    | verbose            MORE TEXT DUDE !!!!!!
    -h    | help               You already know it ...\n\n


Merge all gff3 into one. (Used for permutations test)
    
";
    exit;
}
