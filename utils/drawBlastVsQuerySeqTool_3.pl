#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use Getopt::Long;
use List::Util qw(min max);
use GD::SVG;

###> Loading Options
my %config;
GetOptions (\%config, 'input_file=s', 'type=s', 'id=s', 'ref_length=i', 'round=i', 'sort=s', 'line_width=i', 'overlap', 'verbose', 'help');

#if (
#    $config{help} or 
#    (!exists $config{input_file} or !exists $config{type}) or
#    ($config{type} ne "gff" and $config{type} ne "blast") or
#    ($config{sort} ne "start" and $config{sort} ne "end" and !exists $config{sort})
#    ) {
if (
    $config{help}
    or (!exists $config{input_file} or !exists $config{type})
    or ($config{type} ne "gff" and $config{type} ne "blast")
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
    -overlap         overlap HSPs ? 
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

my %data;
my $nbseq = 0;
my $margeH = 20;
my $margeL = 20;
my $curseurL = $margeL;
my $curseurH = $margeH;
my $epaisseur = $config{'line_width'};

###> Loading Options
print "Loading $config{input_file} ... "        if $config{verbose};

open (FILE, "<$config{'input_file'}");
if ($config{type} eq "blast") {
    while(<FILE>){
        my @tmp = split(/,/, $_);
        $data{$tmp[1]}{'start'} .= $tmp[6]."\t";
        $data{$tmp[1]}{'end'}   .= $tmp[7]."\t";
        $nbseq++;
    }
}
elsif ($config{type} eq "gff"){
    while(<FILE>){
        my @tmp = split(/\t/, $_);
        my @tmp2 = split(/;/, $tmp[8]);
        $tmp2[1] =~ / (\d+) (\d+)/;
        $data{'start'} .= $1."\t";
        $data{'end'}   .= $2."\t";
        $nbseq++;
    }
}
close FILE;
print "OK ! \n" if $config{verbose};
###> Search for picture size
my $largeur = $config{ref_length} + 2 * $margeL;
my $hauteur;
if (!$config{overlap}){$hauteur = 2 * $epaisseur * $nbseq + 2 * $margeH + 2 * 2 * $epaisseur;}
else                  {$hauteur = 2 * $epaisseur * scalar(keys(%data)) + 2 * $margeH + 2 * 2 * $epaisseur;}

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
my @colors = ($red, $green, $forest_green, $blue, $aquamarine, $dark_turquoise, $yellow, $pink, $orange, $brown, $chocolate, $violet, $grey, $black);

# Draw Background
print "Draw background ...\n" if $config{'verbose'};
$image->filledRectangle(0, 0, $largeur, $hauteur, $ligth_grey);



print "Drawing picture ...\n" if $config{verbose};
# Draw RefSeq/Query
$image->filledRectangle($curseurL,                          $curseurH,
                        $curseurL + $config{ref_length},    $curseurH + $epaisseur,
                        $black);
$curseurH = $curseurH + 2 * $epaisseur;

my $index = -1;
my @orderHit;
foreach my $hits (keys(%data)){
    $index++;
    $orderHit[$index]{min}  = min(split(/\t/, $data{$hits}{'start'}));
    $orderHit[$index]{max}  = max(split(/\t/, $data{$hits}{'end'}));
    $orderHit[$index]{name} = $hits;
}

if    ($config{'sort'} eq "end")   {@orderHit = sort{$b -> {'max'} <=> $a -> {'max'} || $b -> {'min'} <=> $a -> {'min'}} @orderHit;}
elsif ($config{'sort'} eq "start") {@orderHit = sort{$a -> {'min'} <=> $b -> {'min'} || $b -> {'max'} <=> $a -> {'max'}} @orderHit;}

my $colorIndex = 0;
foreach my $hits (@orderHit){
    
    print "========> $hits->{name} <========\n";
    print "colorIndex = $colorIndex\n" if $config{verbose};
    
    my @start = split(/\t/, $data{$hits->{name}}{'start'});
    my @end   = split(/\t/, $data{$hits->{name}}{'end'});
    my $round = $config{'round'};
    my @coord;
    
    for(my $i = 0 ; $i <= $#start ; $i++){
        $coord[$i]{'start'} = $start[$i];
        $coord[$i]{'end'}   = $end[$i];
        if    ($config{'sort'} eq "end")   {$coord[$i]{'round'} = int($end[$i]   / $round) * $round;}
        elsif ($config{'sort'} eq "start") {$coord[$i]{'round'} = int($start[$i] / $round) * $round;}
        $coord[$i]{'length'} = $end[$i] - $start[$i];
    }
    
    #Sort datas
    if    ($config{'sort'} eq "end")   {@coord = sort{$b -> {'round'} <=> $a -> {'round'} || $b -> {'length'} <=> $a -> {'length'}} @coord;}
    elsif ($config{'sort'} eq "start") {@coord = sort{$a -> {'round'} <=> $b -> {'round'} || $b -> {'length'} <=> $a -> {'length'}} @coord;}
    
    
    for(my $i = 0 ; $i <= $#coord ; $i++){
        $image->filledRectangle($curseurL + $coord[$i]{'start'},    $curseurH,
                                $curseurL + $coord[$i]{'end'},      $curseurH + $epaisseur,
                                $colors[$colorIndex]);
        $curseurH = $curseurH + 2 * $epaisseur if (!$config{overlap});
        print "\t$curseurL + $start[$i], $curseurH, $curseurL + $end[$i], $curseurH + $epaisseur\n" if $config{verbose};
    }
    $colorIndex++;
    $colorIndex = ($colorIndex == $#colors) ? 0 : $colorIndex;
    
    $curseurH = $curseurH + 2 * $epaisseur if ($config{overlap});
    
}


open(IMG, ">$config{id}.svg");
# make sure we are writing to a binary stream
binmode IMG;
# Convert the image to svg and print it in file
print IMG $image->svg;
close IMG;

print "Finish ! \n" if $config{verbose};
