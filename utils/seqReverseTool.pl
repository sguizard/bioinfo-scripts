#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::SearchIO;
use Data::Dumper;

###> Loading Options
my %config;
GetOptions (\%config, 'input_file=s', 'output_dir=s', 'rm_blast_dir', 'cons_ref=s', 'verbose');

if ($config{help} or !exists $config{input_file} or !exists $config{output_dir}) {
    print
"USAGE : seqReverseTool.pl -input_file inputfile.fa -output_dir directory
    
Options : 
    -i | input_file     multi fasta file (.fa) (Mandatory)
    -o | output_dir     output directory for generated files (Mandatory)
    -cons_ref           file name of the consensus file that you wanna use as reference
    -rm_blast_dir       rm blast dir
    -verbose            MORE text dude !!!!
    -help               you already know it ...\n\n";
    exit;
}

if (!exists $config{input_file} or !exists $config{output_dir}) { die "input_file, output_dir are mandatory ! \n";}
if (! -e $config{input_file}) {die "input_file $config{input_file} not exist ! \n";}
if (! -e $config{output_dir}) {die "output_dir $config{output_dir} not exist ! \n";}

print "input_file : $config{input_file}\n" if $config{verbose};
print "output_dir : $config{output_dir}\n" if $config{verbose};
###> Loading Options

my @sens;
my %antisens;

print "Creating Blast directory ... ";
print "\nmkdir $config{output_dir}/blast\n" if $config{verbose};
`rm -rf $config{output_dir}/blast`;
`mkdir $config{output_dir}/blast`;
print "OK ! \n";

extract_seq(\$config{input_file}, \$config{output_dir}, \@sens, \%antisens, \$config{cons_ref});
run_makeblastdb_blast(\$config{output_dir});
reverse_from_blast(\$config{input_file}, \$config{output_dir}, \@sens, \%antisens);

print "Remove Blast directory ... ";
print "\nrm -rf $config{output_dir}/blast\n" if $config{verbose};
`rm -rf $config{output_dir}/blast` if $config{rm_blast_dir};
print "OK ! \n";


###> Functions
###################################################################################
sub extract_seq {
    my $ref_input_file = shift;
    my $ref_output_dir = shift;
    my $ref_sens = shift;
    my $ref_antisens = shift;
    my $ref_cons_ref = shift;

    print "Extracting sequences ... \n";    
    if (defined($$ref_cons_ref)) {
        my $seqIOin  = Bio::SeqIO->new(-format => 'fasta', -file => "<$$ref_cons_ref");
        my $cons_seq = $seqIOin->next_seq();
        $$ref_sens[0] = $cons_seq->primary_id;
        $$ref_sens[1] = $cons_seq->primary_seq->seq();
        print "cp $$ref_cons_ref $$ref_output_dir/blast/sens.fa\n";
        `cp $$ref_cons_ref $$ref_output_dir/blast/sens.fa`;
        
        $seqIOin  = Bio::SeqIO->new(-format => 'fasta', -file => "<$$ref_input_file");
        while(my $seq = $seqIOin->next_seq()){$$ref_antisens{$seq->primary_id} = $seq->primary_seq->seq();}
        `cp $$ref_input_file $$ref_output_dir/blast/test.fa`;
    }
    else {
        my $seqIOin  = Bio::SeqIO->new(-format => 'fasta', -file => "<$$ref_input_file");
        my $seqIOout = Bio::SeqIO->new(-format => 'fasta', -file => ">$$ref_output_dir/blast/sens.fa");
        
        print "\tSearching longest Seq ... \n";
        my $longest_seq = $seqIOin->next_seq();
        while(my $seq = $seqIOin->next_seq()){
            if ($longest_seq->length < $seq->length) {$longest_seq = $seq;}
        }
        print "\tFound ! \n";
        
        $$ref_sens[0] = $longest_seq->primary_id;
        $$ref_sens[1] = $longest_seq->primary_seq->seq();
        
        $seqIOout->write_seq($longest_seq);
        $seqIOout->DESTROY;
        
        $seqIOout = Bio::SeqIO->new(-format => 'fasta', -file => ">$$ref_output_dir/blast/test.fa");
        $seqIOin  = Bio::SeqIO->new(-format => 'fasta', -file => "<$$ref_input_file");
     
        while(my $seq = $seqIOin->next_seq()){
            next if ($seq->primary_id eq $$ref_sens[0]);
            $$ref_antisens{$seq->primary_id} = $seq->primary_seq->seq();
            $seqIOout->write_seq($seq);
        }
        $seqIOout->DESTROY;
    }
    print "OK !\n";
}
###################################################################################
sub run_makeblastdb_blast{
    my $ref_output_dir = shift;
    print "Running makeblastdb and Blastall ... "; 
    `makeblastdb -in $$ref_output_dir/blast/sens.fa -dbtype nucl`;
    `blastn -query $$ref_output_dir/blast/test.fa -db $$ref_output_dir/blast/sens.fa -out $$ref_output_dir/blast/blast.out -dust no -task blastn -num_threads 8`;
#    `blastall -p blastn -i $$ref_output_dir/blast/test.fa -d $$ref_output_dir/blast/sens.fa -o $$ref_output_dir/blast/blast.out`;
    print "OK !\n"; 
}
###################################################################################
sub reverse_from_blast {
    my $ref_input_dir = shift;
    my $ref_output_dir = shift;
    my $ref_sens = shift;
    my $ref_antisens = shift;

    print "Loading blast results ... ";
    my $blastIO = Bio::SearchIO->new(   -format => 'blast',
                                        -file => "$$ref_output_dir/blast/blast.out");
    print "OK ! \n";
    
    print "Creating output file ... ";
    $$ref_input_dir =~ /([^\/]+)\.fa/;
    my $seqIOout = Bio::SeqIO->new( -format => 'fasta',
                                    -file => ">>$$ref_output_dir/$1_sens_ok.fa");
    #print ">>$$ref_output_dir/$1_sens_ok.fa\n";
    print "OK ! \n"; 

    print "Analysing blast and reversing ... ";
    print "\n" if $config{verbose};
    my $first_seq = Bio::Seq->new(  -seq => $$ref_sens[1],
                                    -id  => $$ref_sens[0]);
    $seqIOout->write_seq($first_seq);
    $first_seq->DESTROY;
    
    my $seqIOerror = Bio::SeqIO->new( -format => 'fasta',
                                      -file => ">>$$ref_output_dir/$1_error.fa");
    
    while (my $result = $blastIO->next_result){
        print "-> Query's Name : ".$result->query_name."\n" if $config{verbose};
        my $hit = $result->next_hit;
        #print Dumper(\$hit); exit;
        if (defined($hit)) {
            
            print "-> Subject's Name : ".$hit->name()."\n" if $config{verbose};
            print "-> Query : \t".$hit->strand('query')."\n" if $config{verbose};
            print "-> Sbject : \t".$hit->strand('sbject')."\n" if $config{verbose};
            
            if ($hit->strand('sbject') == -1){
                my $seq = Bio::Seq->new(-seq => $$ref_antisens{$result->query_name}, -id => $result->query_name);
                $seqIOout->write_seq($seq->revcom());
            }
            else {
                my $seq = Bio::Seq->new(-seq => $$ref_antisens{$result->query_name}, -id => $result->query_name);
                $seqIOout->write_seq($seq);
            }
        }
        else {
            my $seq = Bio::Seq->new(-seq => $$ref_antisens{$result->query_name}, -id => $result->query_name);
            $seqIOerror->write_seq($seq);
        }
    }
    $seqIOout->DESTROY;
    print "OK ! \n"; 
}
###> Functions
###################################################################################














