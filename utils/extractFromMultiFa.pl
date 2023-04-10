#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use Getopt::Long;
use Term::ANSIColor;
use File::Copy;
use Bio::SeqIO;

###> Loading Options
my %config;

GetOptions (\%config, 'in_dir=s', 'out_dir=s', 'method=s', 'number=i', 'from=i', 'to=i', 'help', 'verbose');

if ($config{help} or
   (!exists $config{in_dir} or !exists $config{out_dir} or !exists $config{method} or !exists $config{number})
   ) {
    print color 'red';
    print
"USAGE : extractFromMultiFa.pl -in_dir in/dir -out_dir out/dir -method head -number 300
    
Options : 
    -in_dir         directory containing input files (Mandatory)
    -out_dir        directory for output files (Mandatory)
    -method         head, tail or fromto (Mandatory)
    -number         number of sequence insert in one time (Mandatory if method = head or tail)
    -from           extract from seq num    (Mandatory if method = fromto)
    -to             to seq num              (Mandatory if method = fromto)
    -verbose        MORE text dude !!!!
    -help           you already know it ...\n\n";
    exit;
}

if (!exists $config{in_dir}){die "in_dir is mandatory ! \n";}
if (!exists $config{out_dir}){die "out_dir is mandatory ! \n";}
if (!exists $config{method}){die "method is mandatory ! \n";}
if ($config{method} ne "head" and $config{method} ne "tail" and $config{method} ne "fromto"){die "method must be head OR tail OR fromto \n";}
if ($config{method} eq "head") {
    if (!exists $config{number}){
        die "If method  is head, number must be set ! \n";
    }
}
if ($config{method} eq "tail") {
    if (!exists $config{number}){
        die "If method  is tail, number must be set ! \n";
    }
}
if ($config{method} eq "fromto") {
    if (!exists $config{from} and !exists $config{to}){
        die "If method  is fromto, from and to options must be set\n";
    }
}
if (!-d $config{in_dir}) {die "$config{in_dir} do not exist ! \n";}
if (!-d $config{out_dir}){die "$config{out_dir} do not exist ! \n";}
###> Loading Options


###> Initialize parameters
$| = 1;
###> Initialize parameters

###> Create list of file of $config{in_dir}
opendir(INDIR, $config{in_dir});
my @files = readdir(INDIR);
closedir INDIR;
@files = sort(@files);
shift @files;
shift @files;
###> Create list of file of $config{in_dir}

###> Start message printing
my $count = 0;
print "Processed files : $count / ".($#files + 1);
###> Start message printing

###> Parsing files
foreach (@files){
    
    my $file_size = `grep -c '>' $config{in_dir}/$_`;
    chomp $file_size;
    if ($file_size > $config{number}) {
        
        my $seqIn = Bio::SeqIO->new( -format => 'fasta',
                                     -file => "$config{in_dir}/$_");
        my $outName = $_;
        $outName =~ s/.fa/_$config{number}.fa/;
        my $seqOut = Bio::SeqIO->new(-format => 'fasta',
                                     -file => ">$config{out_dir}/$outName");
        
        if ($config{method} eq "head") {
            for (my $i = 0 ; $i < $config{number} ; $i++){
                my $seq = $seqIn->next_seq;
                $seqOut->write_seq($seq);
            }
        }
        elsif ($config{method} eq "tail") {
            for (my $i = 0 ; $i < $file_size - $config{number} - 1 ; $i++){
                my $seq = $seqIn->next_seq;
            }
            for (my $i = $file_size - $config{number} ; $i < $file_size ; $i++){
                my $seq = $seqIn->next_seq;
                $seqOut->write_seq($seq);
            }
        }
        elsif ($config{method} eq "fromto") {
            for (my $i = 0 ; $i < $file_size - $config{number} - 1 ; $i++){
                my $seq = $seqIn->next_seq;
            }
            for (my $i = $config{from} ; $i <= $config{to} ; $i++){
                my $seq = $seqIn->next_seq;
                $seqOut->write_seq($seq);
            }
        }
    }
    else {
        copy("$config{in_dir}/$_", "$config{out_dir}/$_");
    }
    $count++;
    print "\rProcessed files : $count / ".($#files + 1);
}
print "\n";
###> Parsing files
