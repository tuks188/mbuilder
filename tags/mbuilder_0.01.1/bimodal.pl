#!/usr/local/bin/perl
# bimodal.pl, an ad hoc ellipsoid file constructed for dave saylor's talk...
#
$| = 1;

while ($center = <stdin>) {
  print $center;
  
  @parts = split(" ",$center);
  if ($parts[2] < 0.166667) {
    $small = 1;
    $big = 0;
  }
  if (($parts[2] >= 0.1666667) && ($parts[2] < 0.333333)) {
    $big = ($parts[2] - 0.2)/0.2;
    $small = 1 - $big;
  }
  if (($parts[2] >= 0.333333) && ($parts[2] < 0.6666666)) {
    $big = 1;
    $small = 0;
  }
  if (($parts[2] >= 0.6666666) && ($parts[2] < 0.833333)) {
    $small = ($parts[2] - 0.6)/0.2;
    $big = 1 - $small;
  }
  if ($parts[2] >= 0.833333) {
    $small = 1;
    $big = 0;
  }
  $deviate = rand;
  if ($deviate < $small) {
    print join(" ",&smallAxes()),"\n";
  }
  else {
    print join(" ",&bigAxes()),"\n";
  }
  print "1 0 0\n";
  print "0 1 0\n";
  print "0 0 1\n";
  print rand," ",rand," ",rand,"\n";
}

#
# Small ~20 grain model
#
sub smallAxes {
  return(0.300+rand()*0.050,0.300+rand()*0.050,0.300+rand()*0.050);
}

sub bigAxes {
  return(0.30+rand()*0.05,0.30+rand()*0.05,0.30+rand()*0.05);
}

#
# Medium-small ~150 grain model
#
#sub smallAxes {
#  return(0.100+rand()*0.010,0.100+rand()*0.010,0.100+rand()*0.010);
#}
#
#sub bigAxes {
#  return(0.10+rand()*0.01,0.10+rand()*0.01,0.10+rand()*0.01);
#}
#

#
# Medium  grain model
# 6/7/07  changed to 0.075
#sub smallAxes {
#  return(0.06+rand()*0.005,0.06+rand()*0.005,0.06+rand()*0.005);
#}
#                                                                                
#sub bigAxes {
#  return(0.06+rand()*0.005,0.06+rand()*0.005,0.06+rand()*0.005);
#}

                                                                                




