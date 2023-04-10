#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use List::Util qw(min max);
use GD::SVG;

###> Loading Options
my %config;
GetOptions (\%config, 'input_file=s', 'line_width=i', 'verbose', 'help');

if ($config{help} or !exists $config{input_file}) {
    print
"USAGE : drawMultipleAlignement.pl -input_file multipleAlign.aln -line_width 10
    
Options : 
    -input_file      File to draw (string, blast => outfmt 10)   (Mandatory)
    -line_width      line width                                  (Default: 1)
    -verbose         MORE text dude !!!!
    -help            You already know it ...\n\n";
    exit;
}

if (!exists $config{input_file}) {die "input_file is mandatory ! \n";}
if (!exists $config{line_width}) {$config{line_width} = 1;}
if (! -e    $config{input_file}) {die "input_file $config{input_file} not exist ! \n";}

print "input_file : \t$config{input_file}\n"    if $config{verbose};
print "line_width : \t$config{line_width}\n"    if $config{verbose};

my $epaisseur = $config{'line_width'};
my $length_seq;
my $nbseq = 0;
my $margeH = 20;
my $margeL = 20;
my $curseurL = $margeL;
my $curseurH = $margeH - $epaisseur;



print "Loading $config{input_file} ... "        if $config{verbose};
###> Loading Options

#get length seq
my $seqIOin  = Bio::SeqIO->new(-format => 'fasta', -file => "<$config{'input_file'}");
while (my $seq = $seqIOin->next_seq()){
    $length_seq = $seq->length() if (!$length_seq);
    $nbseq++;
}

print "OK ! \n" if $config{verbose};

###> Search for picture size
my $largeur = $length_seq + 2 * $margeL;
my $hauteur = $epaisseur * $nbseq + 2 * $margeH + 2;


###> Ask if image size is OK
while (1) {
    print "Is the picture size will be ok ? ($largeur px x $hauteur px) [y/n] : ";
    chomp (my $response = <STDIN>);
    if    ($response eq "n") {exit;}
    elsif ($response eq "y") {last;}
}

# Création de l'image
print "Create Picture ...\n" if $config{'verbose'};
my $image = GD::SVG::Image->new($largeur, $hauteur);

# Allocate some colors
my $white =             $image->colorAllocate(255,255,255);
my $red =               $image->colorAllocate(255,0,0);
my $green =             $image->colorAllocate(0,255,0);
my $forest_green =      $image->colorAllocate(34,139,34);
my $blue =              $image->colorAllocate(0,0,255);
my $aquamarine =        $image->colorAllocate(127,255,212);
my $dark_turquoise =    $image->colorAllocate(0,206,209);
my $yellow =            $image->colorAllocate(255,255,0);
my $pink =              $image->colorAllocate(255,192,203);
my $orange =            $image->colorAllocate(255,165,0);
my $brown =             $image->colorAllocate(165,42,42);
my $chocolate =         $image->colorAllocate(210,105,30);
my $violet =            $image->colorAllocate(138,43,226);
my $ligth_grey =        $image->colorAllocate(211,211,211);
my $grey =              $image->colorAllocate(84,84,84);
my $black =             $image->colorAllocate(0,0,0);
#my @colors = ($red, $green, $forest_green, $blue, $aquamarine, $dark_turquoise, $yellow, $pink, $orange, $brown, $chocolate, $violet, $grey, $black);
my %color;
$color{A} = $red;
$color{a} = $red;
$color{C} = $green;
$color{c} = $green;
$color{G} = $yellow;
$color{g} = $yellow;
$color{T} = $blue;
$color{t} = $blue;

# Draw Background
print "Draw background ...\n" if $config{'verbose'};
$image->filledRectangle(0, 0, $largeur, $hauteur, $ligth_grey);



print "Drawing picture ...\n" if $config{verbose};
# Draw 
open(ALN, "<$config{input_file}") or die "Can not open $config{input_file} ! \n";
while (<ALN>) {
    if (/>/){
        $curseurH += $epaisseur;
        $curseurL = $margeL;
        next;
    }
    chomp;
    foreach my $base (split(//, $_)){
        #print "base = $base\n";
        if ($base ne "-"){
            $image->filledRectangle($curseurL,      $curseurH,
                                    $curseurL + 1,  $curseurH + $epaisseur,
                                    $color{$base});
        }
        $curseurL++ ;
    }
}


close ALN;

open(IMG, ">$config{input_file}.svg");
# make sure we are writing to a binary stream
binmode IMG;
# Convert the image to svg and print it in file
print IMG $image->svg;
close IMG;

print "Finish ! \n" if $config{verbose};
