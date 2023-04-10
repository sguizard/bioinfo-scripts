#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use Getopt::Long;
use File::Copy;
use Bio::SeqIO;

###> Loading Options
my %config;

GetOptions (\%config, 'fnodes=s', 'fasta_directory=s', 'verbose', 'help');

if ($config{help} or (!exists $config{fnodes} and !exists $config{fasta_directory})) {
    print
"USAGE : silix-sortSeq.pl -fnodes seq.fnodes -fasta_directory directory_containing_fasta_seq
    
Options : 
    -fnodes             silix .fnodes file (Mandatory)
    -fasta_directory    Directory containing the seq in fasta format
                        fasta's names must be header of sequence without '>' (Mandatory)
    -verbose            MORE text dude !!!!
    -help               you already know it ...\n\n";
    exit;
}

if (!exists $config{fnodes}){die "fnodes is mandatory ! \n";}
if (!exists $config{fasta_directory}){die "fasta_directory is mandatory ! \n";}

if (!-e $config{fnodes}){die "$config{fnodes} do not exist ! \n";}
if (!-d $config{fasta_directory}){die "$config{fasta_directory} do not exist ! \n";}
###> Loading Options

my %cluster;
my %seq;

print "fnodes file :        $config{fnodes}\n"            if $config{verbose};
print "fasta directory :    $config{fasta_directory}\n\n" if $config{verbose};

print "Loading cluster with size ... "      if $config{verbose};
open(FNODE, "<", $config{fnodes}) or die "Can't read $config{fnodes} ! \n";
while (<FNODE>) {
    chomp($_);
    $_ =~ /(.+)\t(.+)/;
    my $cluster_name = $1;
    my $seq_name = lc($2);
    $cluster{$cluster_name}{seq} .= $seq_name."\t";
    $cluster{$cluster_name}{size}++;
    my $inseq = Bio::SeqIO->new(-format => 'fasta', 
                                -file => "<$config{fasta_directory}/$seq_name.fasta");
    $seq{$seq_name} = $inseq->next_seq;
    $inseq->close;
}
close FNODE;
print "OK ! \n" if $config{verbose};

print "Creating cluster directory ... \n"     if $config{verbose};
print "\t-> Creating cluster directory\n"   if $config{verbose};
mkdir("clusters");
print "\t-> Creating not_clustered_seq directory\n"  if $config{verbose};
mkdir("not_clustered_seq");
print "OK ! \n"                             if $config{verbose};

print "Sorting Sequences ... \n"              if $config{verbose};
foreach my $family (keys(%cluster)){
    print "\t-> Processed Family : $family \n" if $config{verbose};
    if ($cluster{$family}{size} > 1) {
        my $outseq = Bio::SeqIO->new(-format => 'fasta', 
                                     -file => ">clusters/$family.fasta");
        foreach (split(/\t/, $cluster{$family}{seq})){
            $outseq->write_seq($seq{$_}) if exists($seq{$_});
        }
        $outseq->close();
    }
    else {
        $cluster{$family}{seq} =~ s/\t//;
        copy("$config{fasta_directory}/$cluster{$family}{seq}.fasta", "not_clustered_seq/$cluster{$family}{seq}.fasta");
    }
}
print "OK ! \n"                             if $config{verbose};

