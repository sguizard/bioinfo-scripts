#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

###> Loading Options
my %config;

GetOptions (\%config, 'input_seq=s', 'sort=s', 'split', 'n=i', 'verbose', 'help');

if ($config{help} or !exists $config{input_seq}) {
    print
"USAGE : sortAlignment.pl -i gal_v6_cons_FAM000003_sens_ok.aln -sort start
    
Options : 
    -i | input_seq  seq file (.fa, .aln) (Mandatory)
    -sort           start (Default) or end or lengthc (croissant) or lengthd (decroissant)
    -split          split alignment
    -n              number of seq by file (Default : 20)
    -verbose        MORE text dude !!!!
    -help           you already know it ...\n\n";
    exit;
}
if (!exists $config{input_seq}){die "input_seq is mandatory ! \n";}
my $seq_by_file;
if (!exists $config{n}){$seq_by_file = 20;}
else {$seq_by_file = $config{n};}
if (!exists $config{sort} or ($config{sort} ne "end" and $config{sort} ne "lengthc" and $config{sort} ne "lengthd")) {$config{sort} = "start";}
print "input file : $config{input_seq}\n"    if $config{verbose};
print "sort by    : $config{sort}\n"         if $config{verbose};
###> Loading Options

my $seqIOin  = Bio::SeqIO->new( -format => 'fasta', 
                                -file => "<$config{input_seq}");

my @aln;
my $index = -1;

while (my $seq = $seqIOin->next_seq()) {
    $index++;
    $aln[$index]{id} = $seq->primary_id;
    $aln[$index]{seq} = $seq->primary_seq->seq();
    $aln[$index]{length} = $seq->primary_seq->length();
    my $seqtmp = $seq->primary_seq->seq();
    $seqtmp =~ s/[actgACGT]+/x/g;
    if ($seqtmp eq "x") {
        $aln[$index]{end} = $seq->primary_seq->length();
        $aln[$index]{start} = 0;
    }
    else {
        my @tmp = split(/x/, $seqtmp);
        $aln[$index]{end} = length($tmp[$#tmp]);
        $aln[$index]{start} = length($tmp[0]);
    }
}

if      ($config{sort} eq "end")      {@aln = sort{$a -> {'end'}    <=> $b -> {'end'}}    @aln;}
elsif   ($config{sort} eq "start")    {@aln = sort{$a -> {'start'}  <=> $b -> {'start'}}  @aln;}
elsif   ($config{sort} eq "lengthc")  {@aln = sort{$a -> {'length'} <=> $b -> {'length'}} @aln;}
elsif   ($config{sort} eq "lengthd")  {@aln = sort{$b -> {'length'} <=> $a -> {'length'}} @aln;}

if ($config{split}) {
    my $count_file = 0; 
    #my $count_seq = 20;
    my $count_seq = $seq_by_file;
    my $seqIOout;
    
    for (my $i = 0 ; $i <= $index ; $i++){
        if ($count_seq == $seq_by_file) {
            $count_seq = 0;
            $count_file++;
            $seqIOout = Bio::SeqIO->new(-format => 'fasta', 
                                        -file => ">$config{input_seq}_sorted_$count_file");
        }
        my $seq = Bio::Seq->new(    -seq => $aln[$i]{seq},
                                    -id =>  $aln[$i]{id});
        $seqIOout->write_seq($seq);
        $count_seq++;
    }
}
else {
    my $seqIOout = Bio::SeqIO->new( -format => 'fasta', 
                                    -file => ">$config{input_seq}_sorted");
    
    for (my $i = 0 ; $i <= $index ; $i++){
        my $seq = Bio::Seq->new(-seq => $aln[$i]{seq},
                                -id =>  $aln[$i]{id});
        $seqIOout->write_seq($seq);
    }
}
