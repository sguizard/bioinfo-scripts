#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use GD;

open(BLAST, "<$ARGV[0]");

# Param :   0 -> name query
#           1 -> size query
#           2 -> name subject
#           3 -> size subject
my @param = split(/\t/, <BLAST>);
chomp($param[3]);

my $l_max;
if ($param[1] > $param[3]){$l_max = $param[1]}
else {$l_max = $param[3]}

my @match;
my $count_match = 0;

for(<BLAST>){
    chomp($_);
    my @line = split(/\t/, $_);
    $match[$count_match][0] = $line[0];
    $match[$count_match][1] = $line[1];
    $match[$count_match][2] = $line[2];
    $match[$count_match][3] = $line[3];
    $count_match++;
}

close BLAST;

my $margeH = 20;
my $margeL = 20;
my $epaisseur = 10;
my $curseurL = $margeL;
my $curseurH = $margeH;

my $largeur = $l_max + 2 * $margeL;
my $hauteur = $epaisseur * ($count_match + 1) * 2 * 2 + 2 * $margeH + 2 * 2 * $epaisseur;
#my $hauteur = 700;

# create a new image
print "Création de l'image $largeur x $hauteur\n";
my $im = new GD::Image($largeur, $hauteur);

# allocate some colors
my @colors;
my $white =             $im->colorAllocate(255,255,255);
my $red =               $im->colorAllocate(255,0,0);
my $green =             $im->colorAllocate(0,255,0);
my $forest_green =      $im->colorAllocate(34,139,34);
my $blue =              $im->colorAllocate(0,0,255);
my $aquamarine =        $im->colorAllocate(127,255,212);
my $dark_turquoise =    $im->colorAllocate(0,206,209);
my $yellow =            $im->colorAllocate(255,255,0);
my $pink =              $im->colorAllocate(255,192,203);
my $orange =            $im->colorAllocate(255,165,0);
my $brown =             $im->colorAllocate(165,42,42);
my $chocolate =         $im->colorAllocate(210,105,30);
my $violet =            $im->colorAllocate(138,43,226);
my $ligth_grey =        $im->colorAllocate(211,211,211);
my $grey =              $im->colorAllocate(84,84,84);
my $black =             $im->colorAllocate(0,0,0);

@colors = ($red, $green, $forest_green, $blue, $aquamarine, $dark_turquoise, $yellow, $pink, $orange, $brown, $chocolate, $violet, $ligth_grey, $grey, $black);

# make the background transparent and interlaced
#$im->transparent($white);
#$im->interlaced('true');

$im->filledRectangle($curseurL, $curseurH,
                     $curseurL + $param[1], $curseurH + $epaisseur,
                     $black);
$curseurH = $curseurH + 2 * $epaisseur;

for (my $i = 0 ; $i < $count_match ; $i++){
    print "$curseurL + $match[$i][0], $curseurH, $curseurL + $match[$i][1], $curseurH + $epaisseur, $colors[$i]\n";
    $im->filledRectangle($curseurL + $match[$i][0], $curseurH,
                         $curseurL + $match[$i][1], $curseurH + $epaisseur,
                         $colors[$i]);
    $curseurH = $curseurH + 2 * $epaisseur;  
}

$curseurH = $curseurH + 2 * 2 * $epaisseur;

$im->filledRectangle($curseurL, $curseurH,
                     $curseurL + $param[3], $curseurH + $epaisseur,
                     $black);
$curseurH = $curseurH + 2 * $epaisseur;

for (my $i = 0 ; $i < $count_match ; $i++){
    $im->filledRectangle($curseurL + $match[$i][2], $curseurH,
                         $curseurL + $match[$i][3], $curseurH + $epaisseur,
                         $colors[$i]);
    $curseurH = $curseurH + 2 * $epaisseur;  
}

open(IMG, ">image.png");
# make sure we are writing to a binary stream
binmode IMG;
# Convert the image to PNG and print it in file
print IMG $im->png;
close IMG;

