#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use Getopt::Long;
#use Term::ANSIColor;

#> Setting Parameters
##> Define outputs colors
#print STDOUT color 'blue';
#print STDERR color 'red';

##> Define options
my %config;
GetOptions (\%config,
            'satInventory=s',
            'verbose',
            'help');

##> Print USAGE if --help
if ($config{help})                      {printUsage(1);}
if (!exists $config{satInventory}) {printError("satInventory option is MANDATORY ! \n", 0); printUsage(1);}

#> Setting parameters
open(SAT, "<$config{satInventory}") or printError("Can not open $config{satInventory} ! ", 1);

my %satHash;
my $count1 = 0;
my $count2 = 0;

while (<SAT>) {
    chomp;
    my ($unit_seq, $unit_size, $ta_num, $repeat_num) = split("\t");
    
    next if ($ta_num < 50 and $unit_size <= 60);
    #next if length($unit_seq) > 10;
    $count1++;
    
    my $reversed_unit = $unit_seq;
    $reversed_unit =~ tr/ACGT/TGCA/;
    $reversed_unit = reverse($reversed_unit);
    
    if (!exists($satHash{$unit_seq}) and !exists($satHash{$reversed_unit})) {
        $satHash{$unit_seq}{ta_num}     = $ta_num;
        $satHash{$unit_seq}{repeat_num} = $repeat_num;
        $satHash{$unit_seq}{unit_size}  = $unit_size;
        $count2++;
    }
    elsif (exists($satHash{$reversed_unit})) {
        $satHash{$reversed_unit}{ta_num}     += $ta_num;
        $satHash{$reversed_unit}{repeat_num} += $repeat_num;
    }
}

open(SIMPLE, ">simpleReport.csv")   or die printError("Can not open simpleReport", 1);
open(MICRO,  ">microsatReport.csv") or die printError("Can not open microsatReport.csv", 1);
open(MINI,   ">minisatReport.csv")  or die printError("Can not open minisatReport", 1);
open(SAT,    ">satReport.csv")      or die printError("Can not open satReport", 1);

print SIMPLE "unit_seq\tunit_size\tta_num\trepeat_num\n";
print MICRO  "unit_seq\tunit_size\tta_num\trepeat_num\n";
print MINI   "unit_seq\tunit_size\tta_num\trepeat_num\n";
print SAT    "unit_seq\tunit_size\tta_num\trepeat_num\n";

foreach my $sat (keys (%satHash)){
    if ($satHash{$sat}{unit_size} == 1) {
        print SIMPLE "$sat\t$satHash{$sat}{unit_size}\t$satHash{$sat}{ta_num}\t$satHash{$sat}{repeat_num}\n";
    }
    elsif ($satHash{$sat}{unit_size} >= 2 and $satHash{$sat}{unit_size} <= 10) {
        print MICRO  "$sat\t$satHash{$sat}{unit_size}\t$satHash{$sat}{ta_num}\t$satHash{$sat}{repeat_num}\n";
    }
    elsif ($satHash{$sat}{unit_size} >= 11 and $satHash{$sat}{unit_size} <= 60) {
        print MINI   "$sat\t$satHash{$sat}{unit_size}\t$satHash{$sat}{ta_num}\t$satHash{$sat}{repeat_num}\n";
    }
    elsif ($satHash{$sat}{unit_size} >= 61) {
        print SAT    "$sat\t$satHash{$sat}{unit_size}\t$satHash{$sat}{ta_num}\t$satHash{$sat}{repeat_num}\n";
    }
    else {printError("You are in the fourth dimension ...\n", 1);}
}

print "total => $count1\tselected =>$count2\n" if $config{verbose};
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
"USAGE : quickCountIsPerWin.pl -gff_t /path/to/gff_te.gff3 -gff_g /path/to/gff_te.gff3 -v
    
Options : 
    -s     | satInventory       file containing the following mysql command output :
                                'use trf;
                                select sequence_unit, unit_size, count(sequence_unit), round(sum(copy_number),0) /
                                from data /
                                group by sequence_unit order by sum(copy_number) desc;'
    -v     | verbose            MORE TEXT DUDE !!!!!!
    -h     | help               You already know it ...\n\n
    
    Remove unit duplicate by eliminating reverse complement sequences
";
    exit;
}