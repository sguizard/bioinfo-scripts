#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use Getopt::Long;
use GD;

###> Loading Options
my %config;
GetOptions (\%config, 'input_file=s', 'type=s', 'id=s', 'ref_length=i', 'round=i', 'sort=s', 'line_width=i', 'verbose', 'help');

#if (
#    $config{help} or 
#    (!exists $config{input_file} or !exists $config{type}) or
#    ($config{type} ne "gff" and $config{type} ne "blast") or
#    ($config{sort} ne "start" and $config{sort} ne "end" and !exists $config{sort})
#    ) {
if (
    $config{help} or 
    (!exists $config{input_file} or !exists $config{type}) or
    ($config{type} ne "gff" and $config{type} ne "blast")
    ) {
    print
"USAGE : drawBlastVsQuerySeqV2.pl -input_file FAM000015_no_match_part.gff3 -type gff -id FAM15 -ref_length 6917 -round 250 -sort start
    
Options : 
    -input_file      File to draw (string, blast => outfmt 10)   (Mandatory)
    -type            Type of the file : blast or gff             (Mandatory)
    -id              Name of the query (string)                  (Mandatory)
    -ref_length      Length of the query (string)                (Mandatory)
    -round           Round factor for sorting (int)              (Default : 500)
    -sort            Sort annotations by : start or end (string) (Default : end)
    -line_width      line width                                  (Default: 10)
    -verbose         MORE text dude !!!!
    -help            You already know it ...\n\n";
    exit;
}

if (!exists $config{input_file}) {die "input_file is mandatory ! \n";}
if (!exists $config{type})       {die "type is mandatory ! \n";}
if (!exists $config{id})         {die "id is mandatory ! \n";}
if (!exists $config{ref_length}) {die "ref_length is mandatory ! \n";}
if (!exists $config{round})      {$config{round} = 500;}
if (!exists $config{sort})       {$config{sort} = "end";}
if (!exists $config{line_width}) {$config{line_width} = 10;}
if (! -e    $config{input_file}) {die "input_file $config{input_file} not exist ! \n";}

print "input_file : \t$config{input_file}\n"    if $config{verbose};
print "type : \t$config{type}\n"                if $config{verbose};
print "id : \t$config{id}\n"                    if $config{verbose};
print "ref_length : \t$config{ref_length}\n"    if $config{verbose};
print "round : \t$config{round}\n"              if $config{verbose};
print "sort : \t\t$config{sort}\n"              if $config{verbose};
###> Loading Options

print "Loading $config{input_file} ... "        if $config{verbose};
my %data;
my $id =     $config{id};
my $length = $config{ref_length};

$data{$id}{'length'} = $length;

open (FILE, "<$config{'input_file'}");
if ($config{type} eq "blast") {
    while(<FILE>){
        my @tmp = split(/,/, $_);
        $data{$id}{'start'} .= $tmp[6]."\t";
        $data{$id}{'end'}   .= $tmp[7]."\t";
        $data{$id}{'nbseq'}++;
    }
}
elsif ($config{type} eq "gff"){
    while(<FILE>){
        my @tmp = split(/\t/, $_);
        my @tmp2 = split(/;/, $tmp[8]);
        $tmp2[1] =~ / (\d+) (\d+)/;
        $data{$id}{'start'} .= $1."\t";
        $data{$id}{'end'}   .= $2."\t";
        $data{$id}{'nbseq'}++;
    }
}
close FILE;
print "OK ! \n" if $config{verbose};

print "Drawing picture ...\n" if $config{verbose};


my @start = split(/\t/, $data{$id}{'start'});
my @end = split(/\t/, $data{$id}{'end'});
my $round = $config{'round'};
my @coord;

for(my $i = 0 ; $i <= $#start ; $i++){
    $coord[$i]{'start'} = $start[$i];
    $coord[$i]{'end'} = $end[$i];
    if    ($config{'sort'} eq "end")   {$coord[$i]{'round'} = int($end[$i]   / $round) * $round;}
    elsif ($config{'sort'} eq "start") {$coord[$i]{'round'} = int($start[$i] / $round) * $round;}
    $coord[$i]{'length'} = $end[$i] - $start[$i];
}

#Sort datas
if    ($config{'sort'} eq "end")   {@coord = sort{$b -> {'round'} <=> $a -> {'round'} || $b -> {'length'} <=> $a -> {'length'}} @coord;}
elsif ($config{'sort'} eq "start") {@coord = sort{$a -> {'round'} <=> $b -> {'round'} || $b -> {'length'} <=> $a -> {'length'}} @coord;}

# definition des variables de l'image
my $margeH = 20;
my $margeL = 20;
#    my $epaisseur = 10;
my $epaisseur = $config{'line_width'};
my $curseurL = $margeL;
my $curseurH = $margeH;

my $largeur = $data{$id}{'length'} + 2 * $margeL;
my $hauteur = 2 * $epaisseur * $data{$id}{'nbseq'} + 2 * $margeH + 2 * 2 * $epaisseur;

# create a new image
#print "Création de l'image $largeur x $hauteur\n";
my $im = new GD::Image($largeur, $hauteur);

# allocate some colors
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
my @colors = ($red, $green, $forest_green, $blue, $aquamarine, $dark_turquoise, $yellow, $pink, $orange, $brown, $chocolate, $violet, $ligth_grey, $grey, $black);

$im->filledRectangle(   $curseurL, $curseurH,
                        $curseurL + $data{$id}{'length'}, $curseurH + $epaisseur,
                        $black);
$curseurH = $curseurH + 2 * $epaisseur;

for(my $i = 0 ; $i <= $#coord ; $i++){
    $im->filledRectangle(   $curseurL + $coord[$i]{'start'}, $curseurH,
                            $curseurL + $coord[$i]{'end'}, $curseurH + $epaisseur,
                            $blue);
    $curseurH = $curseurH + 2 * $epaisseur;
    print "\t$curseurL + $start[$i], $curseurH, $curseurL + $end[$i], $curseurH + $epaisseur\n" if $config{verbose};
}

open(IMG, ">$id.png");
# make sure we are writing to a binary stream
binmode IMG;
# Convert the image to PNG and print it in file
print IMG $im->png;
close IMG;

print "Finish ! \n" if $config{verbose};
