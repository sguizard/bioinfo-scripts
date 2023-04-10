#!/usr/bin/perl -w

use strict;
use diagnostics;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use List::Util qw(min max);
use POSIX;
use Term::ANSIColor;

###> Loading Options
my %config;
GetOptions (\%config,
            'input_file=s',
            'out_file_name=s',
            'line_width=i',
            'marge_top_bottom=i',
            'marge_left_rigth=i',
            'scale_factor=i',
            'scale_factor_auto=i',
            'type_picture=s',
            'auto_yes',
            'verbose',
            'help');
print STDERR color 'red';
$SIG{__WARN__} = sub {};

my $out_file;
my $margeH = 20;
my $margeL = 20;

if    ($config{help} or
       !exists $config{input_file} or
       !exists $config{type_picture})       {printUsage(1);}
if    ($config{type_picture} ne "svg" &&
       $config{type_picture} ne "png")      {printError("Invalid picture type ! \n", 1);}
if    (!exists $config{input_file})         {printError("input_file is mandatory ! \n", 1);}

if    (!exists $config{line_width})         {$config{line_width} = 1;}
if    (!exists $config{scale_factor})       {$config{scale_factor} = 1;}
if    ($config{marge_left_rigth})           {$margeL = $config{marge_left_rigth};}
if    ($config{marge_top_bottom})           {$margeH = $config{marge_top_bottom};}
if    (! -e $config{input_file})            {printError("input_file $config{input_file} not exist ! \n", 1);}
if    ($config{type_picture} eq "svg")      {use GD::SVG;}
elsif ($config{type_picture} eq "png")      {use GD;}
if (!exists $config{out_file_name})         {$out_file = $config{input_file};}
else                                        {$out_file = $config{out_file_name};}
###> Loading Options


print "Input_file : \t$config{input_file}\n" if $config{verbose};
print "Line_width : \t$config{line_width}\n" if $config{verbose};

###> Set variables
my $epaisseur = $config{'line_width'};
my $sf = $config{scale_factor};
my $length_seq;
my $nbseq = 0;
my $curseurL = $margeL;
my $curseurH = $margeH;
my $minus_quantity = 0;
my @density_vector;
my $index_dv = 0;

print "Loading $config{input_file} ... "        if $config{verbose};

###> Get length seq
my $seqIOin  = Bio::SeqIO->new(-format => 'fasta', -file => "<$config{'input_file'}");
while (my $seq = $seqIOin->next_seq()){
    $length_seq = $seq->length() if (!$length_seq);
    $nbseq++;
}
print "OK ! \n" if $config{verbose};
print "Length Seq = $length_seq\n"  if $config{verbose};
print "Number seq = $nbseq\n"       if $config{verbose};


###> Search for picture size
my $hauteur =   $margeH                 # marge haut
              + ($epaisseur * $nbseq)   # épaisseur alignement
              + 10                      # espace entre alignement et densité
              + $epaisseur              # épaisseur density
              + 10                      # espace entre densité et info
              + gdTinyFont->height      # épaisseur info
              + gdTinyFont->height      # légende
              + $margeH;                # marge bas
my $largeur;
# Search for the good scale factor
if (exists $config{scale_factor_auto}) {
    $sf = 1;
    while (1) {
        $largeur = ceil($length_seq/$sf) + 2 * $margeL;
        last if ($largeur < $config{scale_factor_auto});
        $sf++;
    }
}
# Set largeur
else {$largeur = ceil($length_seq/$sf) + 2 * $margeL;}
$largeur = 461 + 2 * $margeL if ($largeur <= 461); # largeur légende


###> Ask if image size is OK
if (!exists $config{auto_yes}){
    while (1) {
        print "Is the picture size will be ok ? ($largeur px x $hauteur px) [y/n] : ";
        chomp (my $response = <STDIN>);
        if    ($response eq "n") {exit;}
        elsif ($response eq "y") {last;}
    }
}


###> Création de l'image
print "Create Picture ...\n" if $config{'verbose'};
my $image;
if      ($config{type_picture} eq "svg") {$image = GD::SVG::Image->new($largeur, $hauteur);}
elsif   ($config{type_picture} eq "png") {$image = GD::Image->new($largeur, $hauteur);}


###> Allocate some colors
#my $red =               $image->colorAllocate(255,0,0);
#my $green =             $image->colorAllocate(0,255,0);
#my $forest_green =      $image->colorAllocate(34,139,34);
#my $blue =              $image->colorAllocate(0,0,255);
#my $aquamarine =        $image->colorAllocate(127,255,212);
#my $dark_turquoise =    $image->colorAllocate(0,206,209);
#my $yellow =            $image->colorAllocate(255,255,0);
#my $pink =              $image->colorAllocate(255,192,203);
#my $orange =            $image->colorAllocate(255,165,0);
#my $brown =             $image->colorAllocate(165,42,42);
#my $chocolate =         $image->colorAllocate(210,105,30);
#my $violet =            $image->colorAllocate(138,43,226);
#my $grey =              $image->colorAllocate(84,84,84);
#my $level4 =            $image->colorAllocate(237,248,251);
#my $level3 =            $image->colorAllocate(178,226,226);
#my $level2 =            $image->colorAllocate(102,194,164);
#my $level1 =            $image->colorAllocate(35,139,69);
my $white =             $image->colorAllocate(255,255,255);
my $black =             $image->colorAllocate(0,0,0);
my $ligth_grey =        $image->colorAllocate(211,211,211);
my $level1 =            $image->colorAllocate(254,240,217);
my $level2 =            $image->colorAllocate(253,204,138);
my $level3 =            $image->colorAllocate(252,141,89);
my $level4 =            $image->colorAllocate(215,48,31);
my $level_density10 =   $image->colorAllocate(165,0,38);
my $level_density9  =   $image->colorAllocate(215,48,39);
my $level_density8  =   $image->colorAllocate(244,109,67);
my $level_density7  =   $image->colorAllocate(253,174,97);
my $level_density6  =   $image->colorAllocate(254,224,139);
my $level_density5  =   $image->colorAllocate(217,239,139);
my $level_density4  =   $image->colorAllocate(166,217,106);
my $level_density3  =   $image->colorAllocate(102,189,99);
my $level_density2  =   $image->colorAllocate(26,152,80);
my $level_density1  =   $image->colorAllocate(0,104,55);
my %color;
$color{l0}    = $white;
$color{l1}    = $level1;
$color{l2}    = $level2;
$color{l3}    = $level3;
$color{l4}    = $level4;
$color{d_l10} = $level_density10;
$color{d_l9}  = $level_density9;
$color{d_l8}  = $level_density8;
$color{d_l7}  = $level_density7;
$color{d_l6}  = $level_density6;
$color{d_l5}  = $level_density5;
$color{d_l4}  = $level_density4;
$color{d_l3}  = $level_density3;
$color{d_l2}  = $level_density2;
$color{d_l1}  = $level_density1;
$color{white} = $white;
$color{black} = $black;


###> Draw Background
print "Draw background ...\n" if $config{verbose};
$image->filledRectangle(0, 0, $largeur, $hauteur, $ligth_grey);


###> Draw 
print "Drawing picture ...\n" if $config{verbose};
open(ALN, "<$config{input_file}") or die "Can not open $config{input_file} ! \n";
my $seq = "";
while (<ALN>) {
    if (/>/ && $seq ne ""){ # When all seq is concatenate and current line is the next header
        $index_dv = 0;      # Initiate index of the density vector
        
        while ($seq =~ s/^([ACGTNacgtn\-]{1,$sf})//) {
            my $subSeq = $1; #extract the number of base corresponding to scale factor
            print "subSeq = $subSeq\n"          if $config{verbose};
            
            my $countMinus = () = $subSeq =~ /\-/g; # count the number of minus in the subseq
            $minus_quantity += $countMinus;
            
            my $level = (floor(($countMinus/$sf)*100) <= 25) ? "l4" : # set color corresponding 
                        (floor(($countMinus/$sf)*100) <= 50) ? "l3" : # to the quatity of minus
                        (floor(($countMinus/$sf)*100) <= 75) ? "l2" :
                        (floor(($countMinus/$sf)*100) <= 99) ? "l1" :"white";
            
            print "countMinus = $countMinus\n"  if $config{verbose};
            print "level = $level\n"            if $config{verbose};
            
            # Draw the pixel
            $image->filledRectangle($curseurL, $curseurH, $curseurL + 1, $curseurH + $epaisseur, $color{$level});
            $curseurL++;
            
            # Compute if the position is filled
            # in other more than 50% for a pixel
            $density_vector[$index_dv]++ if (50 <= ceil((100 - ($countMinus/$sf) * 100)));
            $index_dv++;
        }
        # Positioning cursor on next seq
        $curseurH += $epaisseur;
        $curseurL = $margeL;
        next;
    }
    next if />/; # skip first header
    chomp;       # otherwise remove cariage return
    $seq .= $_;  # and concatenate the seq
}


###> Draw last seq
$index_dv = 0; # same as the upper loop
while ($seq =~ s/^([ACGTNacgtn\-]{1,$sf})//) {
    my $subSeq = $1;
    print "subSeq = $subSeq\n"          if $config{verbose};

    my $countMinus = () = $subSeq =~ /\-/g;
    $minus_quantity += $countMinus;

    my $level = (floor(($countMinus/$sf)*100) <= 25) ? "l4" :
                (floor(($countMinus/$sf)*100) <= 50) ? "l3" :
                (floor(($countMinus/$sf)*100) <= 75) ? "l2" :
                (floor(($countMinus/$sf)*100) <= 99) ? "l1" :"white";

    print "countMinus = $countMinus\n"  if $config{verbose};
    print "level = $level\n"            if $config{verbose};
    
    $image->filledRectangle($curseurL, $curseurH, $curseurL + 1, $curseurH + $epaisseur, $color{$level});
    $curseurL++;
    
    $density_vector[$index_dv]++ if (50 < ceil(100 - (($countMinus/$sf) * 100)));
    $index_dv++;
}
close ALN;


###> Draw density map
$curseurH += $epaisseur + 10; # Positioning cursor for drawing density map ******* epaisseur => se positionne à la fin de l'alignement; 10 espace fixe entre alignment et carte de densité
$curseurL = $margeL;
foreach my $density (@density_vector){
    #print "density = $density\n";
    my $level = (floor(($density/$nbseq)*100) <= 10) ? "d_l1" : # set color corresponding to base density
                (floor(($density/$nbseq)*100) <= 20) ? "d_l2" :
                (floor(($density/$nbseq)*100) <= 30) ? "d_l3" :
                (floor(($density/$nbseq)*100) <= 40) ? "d_l4" :
                (floor(($density/$nbseq)*100) <= 50) ? "d_l5" :
                (floor(($density/$nbseq)*100) <= 60) ? "d_l6" :
                (floor(($density/$nbseq)*100) <= 70) ? "d_l7" :
                (floor(($density/$nbseq)*100) <= 80) ? "d_l8" :
                (floor(($density/$nbseq)*100) <= 90) ? "d_l9" :"d_l10";
    #print "floor(($density/$nbseq)*100) = ".floor(($density/$nbseq)*100)." = $level\n";
    $image->filledRectangle($curseurL, $curseurH, $curseurL + 1, $curseurH + $epaisseur, $color{$level}); # Draw pixel
    $curseurL++;
}


###> Print info
my $filling_rate    = sprintf("%.2f", 100 - (($minus_quantity / ($nbseq * $length_seq)) * 100)); # compute the filling rate of multiple alignement
my $threeSeqDensity = sprintf("%.2f", (3/$nbseq) * 100);                                         # compute dentity percentage for seq for a giving position
$curseurH += $epaisseur + 10; # Positioning cursor for print ing info ******* epaisseur => se positionne à la fin de carte de densité; 10 espace fixe entre carte de densité et infos
$curseurL = $margeL;
$image->string(gdTinyFont, $curseurL, $curseurH, "Lgth seq: ${length_seq}b -- Nb seq: $nbseq -- 1px = ${sf}bp -- Filling % : $filling_rate % -- 3 seq density = $threeSeqDensity % ",$color{'black'});
$curseurH += gdTinyFont->height;


###> Draw legends
$curseurL  = $margeL;
# Legend Multiple Alignement
$image->filledRectangle($curseurL, $curseurH, $curseurL + (gdTinyFont->width * 6), $curseurH + gdTinyFont->height, $color{"l0"});
$image->string(gdTinyFont, $curseurL, $curseurH, " =0% ", $color{'black'});
$curseurL = $curseurL + (gdTinyFont->width * 6);

my $percent = 25;
foreach my $i (1..4){
    $image->filledRectangle($curseurL, $curseurH, $curseurL + (gdTinyFont->width * 6), $curseurH + gdTinyFont->height, $color{"l$i"});
    $image->string(gdTinyFont, $curseurL, $curseurH, "<=$percent% ", $color{'black'});
    $curseurL = $curseurL + (gdTinyFont->width * 6);
    $percent += 25;
}

# Density legende 
$curseurL += 10;
$percent   = 10;

foreach my $i (1..9){
    $image->filledRectangle($curseurL, $curseurH, $curseurL + (gdTinyFont->width * 6), $curseurH + gdTinyFont->height, $color{"d_l$i"});
    $image->string(gdTinyFont, $curseurL, $curseurH, "<=$percent% ", $color{'black'});
    $curseurL = $curseurL + (gdTinyFont->width * 6);
    $percent += 10;
}
$image->filledRectangle($curseurL, $curseurH, $curseurL + (gdTinyFont->width * 6), $curseurH + gdTinyFont->height, $color{"d_l10"});
$image->string(gdTinyFont, $curseurL, $curseurH, "<=$percent%", $color{'black'});


###> Create image
if      ($config{type_picture} eq "svg") {open(IMG, ">$out_file.svg");}
elsif   ($config{type_picture} eq "png") {open(IMG, ">$out_file.png");}
binmode IMG;
if      ($config{type_picture} eq "svg") {print IMG $image->svg;}
elsif   ($config{type_picture} eq "png") {print IMG $image->png;}
close IMG;

print "Finish ! \n" if $config{verbose};


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
    my $exit = shift;
    
    print STDOUT
"USAGE : drawMultipleAlignement.pl -input_file multipleAlign.aln -line_width 10
    
Options : 
    -input_file           Multiple Alignment in fasta format (string)       (Mandatory)
    -out_file_name        Output file name (string)
    -type_picture         svg or png                                        (Mandatory)
    -scale_factor         Set scale factor for large alignment (integer)    (Default: 1)
    -scale_factor_auto    Max picture size, will automatically
                          search for a fitting scale (integer)
    -line_width           Line width (integer)                              (Default: 1)
    -marge_left_rigth     Left and rigth margin in pixel (integer)          (Default: 20)
    -marge_top_bottom     Top and bottom margin in pixel (integer)          (Default: 20)
    -auto_yes             Automatically accept picture size
    -verbose              MORE text dude !!!!
    -help                 You already know it ...

Create a picture of the alignment.

";
    exit if $exit;
}