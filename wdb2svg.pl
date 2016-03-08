#!/usr/bin/perl

# wdb2svg.pl -- convert vector information formatted in plain-text
# world databank format into an SVG file. Details of this format can
# be found at http://www.evl.uic.edu/pape/data/WDB/

# Four projections are currently available:
#   equirectangular -- direct mapping from degrees to pixels
#   mercator -- partial/clipped map, preserves shapes at the expense of size
#   winkel-tripel -- an attempt to compromise on both size and shape
#   orthographic -- hemispherical maps using circular coordinates

# Author: David Eccles (gringer), 2007 <programming@gringer.org>

use strict;
use warnings;
use Math::Trig;
use Math::Trig 'great_circle_distance';
use POSIX;

sub haversin {
    my $angle = shift;
    return((1 - cos($angle))/2);
}

sub sinc {
    my $angle = shift;
    if($angle == 0){ # substitute discontinuity with limit
        return(1);
    } else {
        return(sin($angle)/$angle);
    }
}

sub reversepath {
    my $path = shift;
    return(join("L",reverse(split(/L/,$path))));
}

sub projectpoint {
    my $projection = shift;
    my $projectedx = shift;
    my $projectedy = shift;
    my $originx = shift;
    my $originy = shift;

    my $mercatormax = 85; # degrees to clip mercator projection at
    my $mercatorFactor = 180 / (log(tan(pi/4 + deg2rad($mercatormax)/2)));
    # scaling factor for y component of mercator projection

    my $wtphi0 = acos(2/pi);
    my $wtFactor = 180; # scaling factor for winkel tripel projection
    if($projection ne "equirectangular"){
        if($projection eq "orthographic"){
            my $phi = deg2rad($projectedy - 90);
            my $lambda = deg2rad($projectedx - 180);
            my $phi0 = deg2rad($originy);
            my $lambda0 = deg2rad($originx + 180);
            my $phi1 = deg2rad(-$originy);  # (lambda1,phi1) is directly
            my $lambda1 = deg2rad($originx);    # opposite (lambda0,phi0)
            my $newphi = $phi;
            my $newlambda = $lambda;
            my $top = 0; # false
            # pi/2 adjusts for different origins, see http://perldoc.perl.org/Math/Trig.html
            my $dist0 = great_circle_distance($lambda0, pi/2 - $phi0, $lambda, pi/2 - $phi);
            my $dist1 = great_circle_distance($lambda1, pi/2 - $phi1, $lambda, pi/2 - $phi);
            if($dist0 > pi/2){
                $top = 1; # true
            }
            my $tx = 90 * (cos($phi) * sin($lambda - $lambda0));
            my $ty = 90 * (cos($phi0) * sin($phi) -
                           sin($phi0) * cos($phi) * cos($lambda-$lambda0));
            if($top){
                # mirrored, because the projection is for the inside of the sphere
                $projectedx = 270 - $tx;
            } else {
                $projectedx = $tx;
            }
            $projectedy = $ty + 90;
        }
        if($projection eq "mercator"){
            # two fmods are needed to ensure x > 0
            $projectedx = fmod(fmod($projectedx - $originx, 360) + 360, 360);
            # inverted y because SVG origin is bottom left
            $projectedy = fmod(fmod(90 - $projectedy, 180) + 180, 180);
            if(($projectedy - 90) > $mercatormax){
                $projectedy = $mercatormax + 90;
            }
            if(($projectedy - 90) < (-1 * $mercatormax)){
                $projectedy = -1 * $mercatormax + 90;
            }
            # note: no need to add 90 degrees because this is already 0..180
            $projectedy = log(tan(deg2rad($projectedy)/2)) * $mercatorFactor + 90;
        }
        if($projection eq "winkel-tripel"){
            # two fmods are needed to ensure x > 0
            $projectedx = fmod(fmod($projectedx - $originx, 360) + 360, 360);
            # inverted y because SVG origin is bottom left
            $projectedy = fmod(fmod(90 - $projectedy, 180) + 180, 180);
            my $phi = deg2rad($projectedy - 90);
            my $lambda = deg2rad($projectedx - 180);
            my $alpha = acos(cos($phi)*
                             cos($lambda/2));
            my $x1 = $lambda;
            my $x2 = (2*cos($phi)*sin($lambda/2))/sinc($alpha);
            my $y1 = $phi;
            my $y2 = sin($phi)/sinc($alpha);
            $projectedx = rad2deg(($x1 + $x2) / 2) + 180;
            $projectedy = rad2deg(($y1 + $y2) / 2) + 90;
        }
    } else {
        # two fmods are needed to ensure x > 0
        $projectedx = fmod(fmod($projectedx - $originx, 360) + 360, 360);
        # inverted y because SVG origin is bottom left
        $projectedy = fmod(fmod(90 - $projectedy, 180) + 180, 180);
    }
    return ( $projectedx * 3, $projectedy * 3);
}

sub usage {
  print(STDERR "usage: ./worldmap2svg.pl <input file(s)> [options]\n");
  print(STDERR "\nOther Options:\n");
  print(STDERR "-help                 : Only display this help message\n");
  print(STDERR "-lat <float>          : Set central latitude\n");
  print(STDERR "-long <float>         : Set central longitude\n");
  print(STDERR "-res <float>          : Resolution of line nodes, in degrees\n");
  print(STDERR "-mark <float>,<float> : Mark a point on the map (E/W, N/S)\n");
  print(STDERR "-projection <name>    : Change projection function, (x,y) = f(lat,long)\n");
  print(STDERR "   currently supports 'equirectangular', 'orthographic',\n");
  print(STDERR "  'mercator', 'winkel-tripel' (default)\n");
  print(STDERR "\n");
}



my $resolution = 0.7; # resolution of line nodes, in degrees
my $joinRes = 0.04; # path ends must be at least this close to be merged

my @pathStartx = ();
my @pathStarty = ();
my @pathEndx = ();
my @pathEndy = ();
my @paths = ();
my $draw = 0;
my $pointsToGo = 0;
my $segment = 0;
my $rank = 0;
my $oldx = 0;
my $oldxorig = 0;
my $oldy = 0;
my $oldyorig = 0;

my @markedx = ();
my @markedy = ();


my $centrex = 0; # East-West centre for image
my $centrey = 0;   # North-South centre for image
                   # (not quite sure how well this works yet)

my $degreecalc = 2*pi / 360;

my $projectionName = "winkel-tripel";
# The script currently supports the following projections
# "equirectangular" -- straight degrees to pixels conversion
# "orthographic" -- two hemisphere circular projection
# "mercator" -- shape preserving, higher latitudes scaled more
# "winkel-tripel" -- mean of equirectangular and aitoff projection

my $makeEdge = 1; # create a file for the border around the edge of the map

my $currentPath = "";
my $segnumber = 0;
my $onlyOnePoint = 0;

my @filenames = ();

# extract command line arguments
while(@ARGV){
    my $argument = shift @ARGV;
    if(-f $argument){ # file existence check
            push(@filenames, $argument);
    } else {
        if($argument eq "-help"){
            usage();
            exit(0);
        }
        elsif($argument eq "-mark"){
            my $location = shift @ARGV;
            if($location =~ /([0-9\.\-]+),([0-9\.\-]+)/){
                my $ypos = $1;
                my $xpos = $2;
                push(@markedx, $xpos);
                push(@markedy, $ypos);
                printf(STDERR "Creating marked feature at ($ypos, $xpos)\n");
            } else {
                printf(STDERR "\nError: Feature location cannot be understood, '$location'\n");
                printf(STDERR "Expecting <float>,<float> (e.g. '-41,174')\n");
                usage();
                exit(7);
            }
        }
        elsif($argument eq "-res"){
            $resolution = shift @ARGV;
            printf(STDERR "setting resolution to $resolution degrees\n");
        }
        elsif($argument eq "-lat"){
            $centrey = shift @ARGV;
            printf(STDERR "setting North/South (latitude) centre to $centrey degrees\n");
        }
        elsif($argument eq "-long"){
            $centrex = shift @ARGV;
            printf(STDERR "setting East/West (longitude) centre to $centrex degrees\n");
        }
        elsif($argument eq "-projection"){
            $projectionName = shift @ARGV;
            printf(STDERR "Using '$projectionName' projection\n");
        }
       else{
            printf(STDERR "Error: Unknown command-line option, '$argument'\n");
            usage();
            exit(4);
        }
    }
}

if($projectionName !~ /^(equirectangular|orthographic|mercator|winkel-tripel)$/){
    printf(STDERR "\nError: Unknown projection, '$projectionName'\n\n");
    usage();
    exit(5);
}

if(($centrey != 0) && ($projectionName !~ /^orthographic$/)){
    printf(STDERR "\nError: Altering central latitude is only allowed for\n");
    printf(STDERR "      'orthographic' projection.\n\n");
    usage();
    exit(6);
}


@ARGV = @filenames;

my $rounddecimals = int(-1 * log($resolution * $resolution)/log(10) + 3);
# print "rounding to $rounddecimals decimal places\n";

# adjust for 2D distance, and scaling of final picture
$resolution = ($resolution * $resolution) * 3;
# adjust for scaling of final picture
$joinRes = $joinRes * 3;

my $lineCount = 0;

print(STDERR "Retrieving paths (one '.' per 100,000 lines)");

while(<>){
    $lineCount++;
    if($lineCount % 100000 == 0){
        print (STDERR ".");
    }
    if(/segment ([0-9]+)  rank ([0-9]+)  points ([0-9]+)/){
        $segment = $1; $rank = $2; $pointsToGo = $3;
        $segnumber++;
        if($rank == 1){
            $paths[$segnumber] = "";
            $onlyOnePoint = ($pointsToGo == 1); # special case, one point found
        }
    }
    elsif(/\s+([0-9\-\.]+) ([0-9\-\.]+)/){
        my $newx = $2;
        my $newy = $1;
        # two fmods are needed to ensure x > 0
        $newx = fmod(fmod($newx, 360) + 360, 360);
        $newy = $newy + 90;
        $oldxorig = $newx;
        $oldyorig = $newy;
        ( $newx, $newy ) = projectpoint($projectionName, $newx, $newy,
                                        $centrex, $centrey);
        $newx = sprintf("%.${rounddecimals}f", $newx);
        $newy = sprintf("%.${rounddecimals}f", $newy);
        $pointsToGo--;
        if($rank == 1){
            my $diff = sqrt(abs(($oldx-$newx))+(abs($oldy-$newy)));
            if($pointsToGo || $onlyOnePoint){
                if($diff>$resolution){
                    if(!$currentPath){
                        if(!$pathStartx[$segnumber]){
                            if(!$newx || !$newy){
                                warn "malformatted variable(s): ($newx, $newy)";
                            }
                            $pathStartx[$segnumber] = $newx;
                            $pathStarty[$segnumber] = $newy;
                        }
                        if($onlyOnePoint){
                            $paths[$segnumber] =
                                qq{<circle cx="${newx}" cy="${newy}"};
                        } else {
                            $paths[$segnumber] =  $paths[$segnumber].
                                qq{${newx},${newy}}."L";
                        }
                    }
                    #print qq{${newx},${newy}L\n};
                    $oldx = $newx;
                    $oldy = $newy;
                }
            }
            if(!$pointsToGo || $onlyOnePoint){
                if(!$newx || !$newy){
                    warn "malformatted variable(s): ($newx, $newy)";
                }
                if($onlyOnePoint){
                    $pathEndx[$segnumber] = -9999;
                    $pathEndy[$segnumber] = -9999;
                } else {
                    $pathEndx[$segnumber] = $newx;
                    $pathEndy[$segnumber] = $newy;
                    $paths[$segnumber] =  $paths[$segnumber].
                        qq{${newx},${newy}};
                    #print qq{$newx,$newy" id="$segnumber"/>\n};
                }
                $oldx = -9999;
                $oldxorig = -9999;
                $oldy = -9999;
                $oldyorig = -9999;
            }
        }
    }
}

if($makeEdge && $projectionName ne "orthographic"){
    $segnumber++;
    $paths[$segnumber] = "";
    my $ressqrt = sqrt($resolution);
    my ( $newx, $newy ) = projectpoint($projectionName, 0, 0, 0, 0);
    $newx = sprintf("%.${rounddecimals}f", $newx);
    $newy = sprintf("%.${rounddecimals}f", $newy);
    $paths[$segnumber] =  $paths[$segnumber].
        qq{${newx},${newy}}."L";
    $pathStartx[$segnumber] = $newx;
    $pathStarty[$segnumber] = $newy;
    for(my $edgex = 0; $edgex < 360; $edgex++){
        ( $newx, $newy ) = projectpoint($projectionName, $edgex, 0, 0, 0);
        $newx = sprintf("%.${rounddecimals}f", $newx);
        $newy = sprintf("%.${rounddecimals}f", $newy);
        $paths[$segnumber] =  $paths[$segnumber].
            qq{${newx},${newy}}."L";
    }
    for(my $edgey = 0; $edgey < 180; $edgey++){
        ( $newx, $newy ) = projectpoint($projectionName, 360-$resolution, $edgey, 0, 0);
        $newx = sprintf("%.${rounddecimals}f", $newx);
        $newy = sprintf("%.${rounddecimals}f", $newy);
        $paths[$segnumber] =  $paths[$segnumber].
            qq{${newx},${newy}}."L";
    }
    for(my $edgex = 360; $edgex > 0; $edgex--){
        ( $newx, $newy ) = projectpoint($projectionName, $edgex, 180-$resolution, 0, 0);
        $newx = sprintf("%.${rounddecimals}f", $newx);
        $newy = sprintf("%.${rounddecimals}f", $newy);
        $paths[$segnumber] =  $paths[$segnumber].
            qq{${newx},${newy}}."L";
    }
    for(my $edgey = 180; $edgey > 0; $edgey--){
        ( $newx, $newy ) = projectpoint($projectionName, 0, $edgey, 0, 0);
        $newx = sprintf("%.${rounddecimals}f", $newx);
        $newy = sprintf("%.${rounddecimals}f", $newy);
        $paths[$segnumber] =  $paths[$segnumber].
            qq{${newx},${newy}}."L";
    }
    ( $newx, $newy ) = projectpoint($projectionName, 0, 0, 0, 0);
    $newx = sprintf("%.${rounddecimals}f", $newx);
    $newy = sprintf("%.${rounddecimals}f", $newy);
    $paths[$segnumber] =  $paths[$segnumber].
        qq{${newx},${newy}};
    $pathEndx[$segnumber] = $newx;
    $pathEndy[$segnumber] = $newy;
}

print(STDERR " done!\n");

print(STDERR "Merging paths...");


# loop through and find path ends that are close enough to each other to
# be merged

my $changed; # method is recursive this
do{
    $changed = 0;
    for(my $i = 0; $i < (@paths + 0); $i++){
        if(($paths[$i]) && (($pathStartx[$i] ne $pathEndx[$i])
                            || ($pathStarty[$i] ne $pathEndy[$i]))){
            if(!$pathEndx[$i] || !$pathEndy[$i] ||
               !$pathStartx[$i] || !$pathStarty[$i]){
#		warn "uninitialised value!!! (position $i):".
#		    "($pathEndx[$i],$pathEndy[$i]),".
#		    "($pathStartx[$i],$pathStarty[$i])";
                next;
            }
            if($pathEndx[$i] != -9999){ # check to make sure it's not a single point
# 	    warn "checking ($pathStartx[$i],$pathStarty[$i]) and ".
# 	         "($pathEndx[$i], and $pathEndy[$i])";
                for(my $o = $i+1; $o < (@paths + 0); $o++){
                    if($paths[$o]){
                        if(!$pathEndx[$o] || !$pathEndy[$o] ||
                           !$pathStartx[$o] || !$pathStarty[$o]){
#			warn "uninitialised value!!! (position $o):".
#			    "($pathEndx[$o],$pathEndy[$o]),".
#			    "($pathStartx[$o],$pathStarty[$o])";
                            next;
                        }
                        #warn "diff (".abs($pathStartx[$o]-$pathEndx[$i]).",".
                        #	abs($pathStarty[$o]-$pathEndy[$i]).")".
                        #	" (".abs($pathStartx[$i]-$pathEndx[$o]).",".
                        #	abs($pathStartx[$i]-$pathEndx[$o]).")";
                        if((abs($pathEndx[$i] - $pathStartx[$o])<$joinRes) &&
                           (abs($pathEndy[$i] - $pathStarty[$o])<$joinRes)){
                            $changed = 1;
                            $paths[$i] = $paths[$i]."L\n".$paths[$o];
                            $paths[$o] = "";
                            $pathEndx[$i] = $pathEndx[$o];
                            $pathEndy[$i] = $pathEndy[$o];
                            $pathStartx[$o] = "";
                            $pathStarty[$o] = "";
                            $pathEndx[$o] = "";
                            $pathEndy[$o] = "";
# 			warn "created new path ".
# 			    "($pathStartx[$i],$pathStarty[$i]) to ".
# 			    "($pathEndx[$i],$pathEndy[$i])";
                        }
                        elsif((abs($pathEndx[$o] - $pathStartx[$i])<$joinRes) &&
                              (abs($pathEndy[$o] - $pathStarty[$i])<$joinRes)){
#			warn "found similar end to path $i";
                            $changed = 1;
                            $paths[$i] = $paths[$o]."L\n".$paths[$i];
                            $paths[$o] = "";
                            $pathStartx[$i] = $pathStartx[$o];
                            $pathStarty[$i] = $pathStarty[$o];
                            $pathStartx[$o] = "";
                            $pathStarty[$o] = "";
                            $pathEndx[$o] = "";
                            $pathEndy[$o] = "";
# 			warn "created new path ".
# 			    "($pathStartx[$i],$pathStarty[$i]) to ".
# 			    "($pathEndx[$i],$pathEndy[$i])";
                        }
                        elsif((abs($pathEndx[$i] - $pathEndx[$o])<$joinRes) &&
                              (abs($pathEndy[$i] - $pathEndy[$o])<$joinRes)){
                            $changed = 1;
                            $paths[$i] = $paths[$i]."L".reversepath($paths[$o]);
                            $paths[$o] = "";
                            $pathEndx[$i] = $pathStartx[$o];
                            $pathEndy[$i] = $pathStarty[$o];
                            $pathStartx[$o] = "";
                            $pathStarty[$o] = "";
                            $pathEndx[$o] = "";
                            $pathEndy[$o] = "";
# 			warn "created new path ".
# 			    "($pathStartx[$i],$pathStarty[$i]) to ".
# 			    "($pathEndx[$i],$pathEndy[$i])";
                        }
                        elsif((abs($pathStartx[$i] - $pathStartx[$o])<$joinRes) &&
                              (abs($pathStarty[$i] - $pathStarty[$o])<$joinRes)){
                            $changed = 1;
                            $paths[$i] = reversepath($paths[$o])."L".$paths[$i];
                            $paths[$o] = "";
                            $pathStartx[$i] = $pathEndx[$o];
                            $pathStarty[$i] = $pathEndy[$o];
                            $pathStartx[$o] = "";
                            $pathStarty[$o] = "";
                            $pathEndx[$o] = "";
                            $pathEndy[$o] = "";
# 			warn "created new path ".
# 			    "($pathStartx[$i],$pathStarty[$i]) to ".
# 			    "($pathEndx[$i],$pathEndy[$i])";
                        }
                    }
                }
            }
        }
    }
} while($changed);

# close paths where start/end is the same
for(my $i = 0; $i < (@paths + 0); $i++){
    if($paths[$i]){
        if((abs($pathEndx[$i] - $pathStartx[$i])<$joinRes) &&
           (abs($pathEndy[$i] - $pathStarty[$i])<$joinRes)){
            $paths[$i] = $paths[$i]."Z";
            if($paths[$i] =~ /^[^L]*L[^L]*$/){ # only one point
                $paths[$i] =~ s/L.*$//;
                $paths[$i] =~ s/^/<circle cx="/;
                $paths[$i] =~ s/,/" cy="/;
                $paths[$i] =~ s/$/"/;
            }
        }
    }
}

print(STDERR " done!\n");

print(STDERR "Generating SVG file...");

# and now the printing begins...

print qq{<?xml version="1.0" standalone="no"?>\n};
print qq{<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"\n};
print qq{"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n};
print qq{<svg width="1080" height="540"\n};
print qq{xmlns="http://www.w3.org/2000/svg" version="1.1">\n};
print qq{<title>World Databank II 2 SVG converter</title>\n};
print qq{<desc>World Map Image derived from WDB, };
print qq{$projectionName projection, using David Eccles };
print qq{(gringer)'s converter</desc>\n};

my $outnumber = 0;
for(@paths){
    if($_){
        $outnumber++;
        if(m/circle/){
            print $_;
            print qq{ r="0.5" style="stroke:none;fill:#000000}
        } else {
            print qq{<path style="fill:none;stroke:#000000;};
            print qq{stroke-width:0.5;stroke-linecap:round;};
            print qq{stroke-linejoin:round" d="M\n};
            s/L/L\n/g;
            print $_;
        }
        print qq{" id="$outnumber"/>\n};
    }
}

if(($projectionName eq "orthographic") && ($makeEdge)){
    print qq{<circle cx="270" cy="270" r="270" style="stroke-width:1;stroke:#000000;fill:none;"};
    print qq{ id="frameFront"/>};
    print qq{<circle cx="810" cy="270" r="270" style="stroke-width:1;stroke:#000000;fill:none;"};
    print qq{ id="frameFront"/>};
}

if(scalar(@markedx)){
    for(my $i=0; $i < scalar(@markedx); $i++){
        my $xpos = fmod(fmod($markedx[$i], 360) + 360, 360);
        my $ypos = $markedy[$i] + 90;
        my ( $xpos, $ypos ) =
            projectpoint($projectionName, $xpos, $ypos, $centrex, $centrey);
        print qq{<circle cx="$xpos" cy="$ypos" r="4" style="stroke-width:1;stroke:#FF0000;fill:none;"};
        print qq{ id="MP$i"/>};
    }
}

print qq{</svg>\n};

print(STDERR " done!\n");
