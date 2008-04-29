#!/usr/local/bin/perl

use IPC::Open2;
#use Getopt::Long;

$| = 1;

if ($#ARGV < 0) {
  print <<END_OF_USAGE;

  Usage::recursiveSampler.pl [-fraction frac] [-bbox xmin xmax ymin ymax zmin zmax] 
                             ellipsoidFunction

    where:
          -fraction frac specifies the threshold of the shortest ellipsoid dimension 
	    below which refinement is triggered.  default is 0.333

	  -bbox specifies the bounding box of the region to be sampled.  0 1 0 1 0 1 by default.

	  ellipsoidFunction points to an executable that reads triples from stdin and writes
	    ellipsoids to stdout.  ellipsoids are defined in the following manner:

	       center (3 doubles)
	       semi-axis lengths (3 doubles)
	       major axis (3 doubles - a unit vector)
               middle axis (3 doubles - a unit vector)
	       minor axis (3 doubles - a unit vector)
	       orientation (3 doubles)

 recursiveSampler samples the computational space (either specified by -bbox, or the unit cube) 
 randomly, and calls ellipsoidFunction at each of the sample points.  it keeps track of the smallest
 dimension of the ellipsoid, and compares that with the dimension of it\'s bounding box.  if the 
 bounding box is greater than the -fraction parameter (0.333 by default) times the smallest dimension
 of the ellipsoid, then recursiveSampler resamples the box with recursively subdivided boxes.
 recursiveSampler prints all the ellipsoids on stdout...

END_OF_USAGE
  exit(1);
}


$fraction = 1.05;
$bbox = [0,1,0,1,0,1];
$GLOBAL_MAX = 1000000;
$listFile = "ellipsoids.txt";
$listbool = 0;
$scale=1;

# foreach $ind (0..$#ARGV) {
# 	print $ARGV[$ind];
# 	print "\n";
# }

# $result = GetOptions ("age=i" => $age);  
# $result = GetOptions("fraction=i" => $fraction);
# print $fraction;
# print "\n";

foreach $ind (0..$#ARGV) {
	
	if ($ARGV[$ind] =~ /^-fraction/) {
		$fraction = $ARGV[$ind+1];
	}
	if ($ARGV[$ind] =~ /^-bbox/) {
		$bbox = [@ARGV[$ind+1..$ind+6]];
	}
	if ($ARGV[$ind] =~ /^-list/) {
		$listbool = 1;
		$listFile = $ARGV[$ind+1];
	} 
}


if ($listbool == 1 ) {
# Open the data file and extract the data set
	open(DAT, $listFile) || die("Could not open file!");
	@raw_data=<DAT>;
	close(DAT);
	$globalCount = 0;
# Recursively sample the list and put ellipsoids in the unit RVE
	&sampleList($bbox,$fraction,$scale);
	
} else {
	$ellipsoidFunction = pop(@ARGV);
	
	# ok, open up a double ended pipe to ellipsoidFunction...
 	$pid = open2(*rdfh,*wrfh,"$ellipsoidFunction");
 	if ($! =~ /^open2/) {
 		die $!;
 	}
	
	
 	&sample($bbox,$fraction,rdfh,wrfh);
	
	close($rdfh);
	close($wrfh);
}

sub sample {
  my $bbox = shift(@_);
  my $fraction = shift(@_);
  my $rdfh = shift(@_);
  my $wrfh = shift(@_);

  my $x;
  my $y;
  my $z;
  my $maxDim;
  my @parts;
    

  $x = $bbox->[0]+(rand())*($bbox->[1] - $bbox->[0]);
  $maxDim = $bbox->[1] - $bbox->[0];

  $y = $bbox->[2]+(rand())*($bbox->[3] - $bbox->[2]);
  if (($bbox->[3] - $bbox->[2]) > $maxDim) {
    $maxDim = $bbox->[3] - $bbox->[2];
  }

  $z = $bbox->[4]+(rand())*($bbox->[5] - $bbox->[4]);
  if (($bbox->[5] - $bbox->[4]) > $maxDim) {
    $maxDim = $bbox->[5] - $bbox->[4];
  }

  print $wrfh "$x $y $z\n";

  my $center = <$rdfh>;
  my $axes = <$rdfh>;
  my $axis1 = <$rdfh>;
  my $axis2 = <$rdfh>;
  my $axis3 = <$rdfh>;
  my $orientation = <$rdfh>;

  print $center;
  print $axes;
  print $axis1;
  print $axis2;
  print $axis3;
  print $orientation;
  $globalCount++;


    

  @parts = split(" ",$axes);
  
#
# Redefine the limiting dimension as the Middle dimension not the Min.
#

#  $minAxis = $parts[0];
#  if ($parts[1] < $minAxis) {
#    $minAxis = $parts[1];
#  }
#  if ($parts[2] < $minAxis) {
#    $minAxis = $parts[2];
#  }
#  if ($minAxis eq  $parts[1]) {
#    if ( $parts[0] < $parts[2] )
#    {
#      $midAxis = $parts[0];
#    } else {
#      $midAxis = $parts[2];
#    }
#  }
#  if ($minAxis eq  $parts[0]) {
#    if ( $parts[1] < $parts[2] )
#    {
#       $midAxis = $parts[1];
#    } else {
#       $midAxis = $parts[2];
#    }
#  }
  @parts = sort(@parts);
  $midAxis = $parts[1];

  if ($maxDim > $fraction*$midAxis) {
    my $midx = ($bbox->[0] + $bbox->[1])/2.0;
    my $midy = ($bbox->[2] + $bbox->[3])/2.0;
    my $midz = ($bbox->[4] + $bbox->[5])/2.0;
  
    &sample([$bbox->[0],$midx,$bbox->[2],$midy,$bbox->[4],$midz],$fraction,$rdfh,$wrfh);
    &sample([$midx,$bbox->[1],$bbox->[2],$midy,$bbox->[4],$midz],$fraction,$rdfh,$wrfh);
    &sample([$bbox->[0],$midx,$midy,$bbox->[3],$bbox->[4],$midz],$fraction,$rdfh,$wrfh);
    &sample([$midx,$bbox->[1],$midy,$bbox->[3],$bbox->[4],$midz],$fraction,$rdfh,$wrfh);
    &sample([$bbox->[0],$midx,$bbox->[2],$midy,$midz,$bbox->[5]],$fraction,$rdfh,$wrfh);
    &sample([$midx,$bbox->[1],$bbox->[2],$midy,$midz,$bbox->[5]],$fraction,$rdfh,$wrfh);
    &sample([$bbox->[0],$midx,$midy,$bbox->[3],$midz,$bbox->[5]],$fraction,$rdfh,$wrfh);
    &sample([$midx,$bbox->[1],$midy,$bbox->[3],$midz,$bbox->[5]],$fraction,$rdfh,$wrfh);
  }

  #if ( $globalCount > $GLOBAL_MAX ){
  #exit(1);
  #}
}


#waitpid $pid, 0;

sub sampleList {
  my $bbox = shift(@_);
  my $fraction = shift(@_);
  my $scale = shift(@_);

  my $x;
  my $y;
  my $z;
  my $maxDim;
  my @parts;


  $x = $bbox->[0]+(rand())*($bbox->[1] - $bbox->[0]);
  $maxDim = $bbox->[1] - $bbox->[0];

  $y = $bbox->[2]+(rand())*($bbox->[3] - $bbox->[2]);
  if (($bbox->[3] - $bbox->[2]) > $maxDim) {
    $maxDim = $bbox->[3] - $bbox->[2];
  }

  $z = $bbox->[4]+(rand())*($bbox->[5] - $bbox->[4]);
  if (($bbox->[5] - $bbox->[4]) > $maxDim) {
    $maxDim = $bbox->[5] - $bbox->[4];
  }

  my @el = randomEllipse();

  my @Axis = split(/\s+/,$el[0]);
  my $a = ("$Axis[0]" / $scale);
  my $b = ("$Axis[1]" / $scale);
  my $c = ("$Axis[2]" / $scale);
  my $center = "$x $y $z\n";
  my $axes = "$a $b $c\n";
  my $axis1 = $el[1];
  my $axis2 = $el[2];
  my $axis3 = $el[3];
  my $orientation = $el[4];

  print $center;
  print $axes;
  print $axis1;
  print $axis2;
  print $axis3;
  print $orientation;
  $globalCount++;

  @parts = split(" ",$axes);
  
#
# Redefine the limiting dimension as the Middle dimension not the Min.
#

#  $minAxis = $parts[0];
#  if ($parts[1] < $minAxis) {
#    $minAxis = $parts[1];
#  }
#  if ($parts[2] < $minAxis) {
#    $minAxis = $parts[2];
#  }
#  if ($minAxis eq  $parts[1]) {
#    if ( $parts[0] < $parts[2] )
#    {
#      $midAxis = $parts[0];
#    } else {
#      $midAxis = $parts[2];
#    }
#  }
#  if ($minAxis eq  $parts[0]) {
#    if ( $parts[1] < $parts[2] )
#    {
#       $midAxis = $parts[1];
#    } else {
#       $midAxis = $parts[2];
#    }
#  }
  @parts = sort(@parts);
  $midAxis = $parts[0];

  if ($maxDim > $fraction*$midAxis) {
    my $midx = ($bbox->[0] + $bbox->[1])/2.0;
    my $midy = ($bbox->[2] + $bbox->[3])/2.0;
    my $midz = ($bbox->[4] + $bbox->[5])/2.0;

    &sampleList([$bbox->[0],$midx,$bbox->[2],$midy,$bbox->[4],$midz],$fraction,$scale);
    &sampleList([$midx,$bbox->[1],$bbox->[2],$midy,$bbox->[4],$midz],$fraction,$scale);
    &sampleList([$bbox->[0],$midx,$midy,$bbox->[3],$bbox->[4],$midz],$fraction,$scale);
    &sampleList([$midx,$bbox->[1],$midy,$bbox->[3],$bbox->[4],$midz],$fraction,$scale);
    &sampleList([$bbox->[0],$midx,$bbox->[2],$midy,$midz,$bbox->[5]],$fraction,$scale);
    &sampleList([$midx,$bbox->[1],$bbox->[2],$midy,$midz,$bbox->[5]],$fraction,$scale);
    &sampleList([$bbox->[0],$midx,$midy,$bbox->[3],$midz,$bbox->[5]],$fraction,$scale);
    &sampleList([$midx,$bbox->[1],$midy,$bbox->[3],$midz,$bbox->[5]],$fraction,$scale);
  }

  #if ( $globalCount > $GLOBAL_MAX ){
  #exit(1);
  #}
}


sub randomEllipse{
my $rline = rand(@raw_data);

# Determine the starting line of the dataset.
my $dline = ($rline / 6);
my $il = int ($dline);

# Read the dataset from the file
my $line1 =  $raw_data[$il*6];
my $line2 =  $raw_data[$il*6+1];
my $line3 =  $raw_data[$il*6+2];
my $line4 =  $raw_data[$il*6+3];
my $line5 =  $raw_data[$il*6+4];
my $line6 =  $raw_data[$il*6+5];

my @lines = ($line2, $line3, $line4, $line5, $line6);

}
    
  
 

