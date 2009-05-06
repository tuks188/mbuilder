#!/usr/local/bin/perl
#  dervied from ....
# bimodal.pl, an ad hoc ellipsoid file constructed for dave saylor's talk...

# intended to provide a log-normal distribution of sizes
#  approximately that observed experimentally for single phase polycrystals
#  ADR, Dec 06

# Updated to read command line size if provided
# if no commandline argument is given then scale = 1;
my $scale = shift;
$scale>0 or $scale=1;

# Lognormal parameters
my $Xspread = 0.6 ;
my $Xoffset = -0.11 ;
my $Yspread = 0.75;
my $Yoffset = -0.11;

# The size of the grain in the Z direction divided by the total size
# of the RVE in the z direction
my $size = 10.0/$scale ; 
#  adjust this parameter to get different sizes of ellipsoids
#  0.1  gave ~8500 ellipsoids  and 151 active cells
#  0.05  gave 71000 ellipsoids and 898 active cells


my $aspect = $size * 5 ;
#  adjust this parameter to vary the range of aspect ratios

$| = 1;

my $i = 0 ;
my $t1 = 0 ;
my $t2 = 0 ;
my $t3 = 0 ;

while ($center = <stdin>) {
  print $center;

#######################################################################
#  scale the size of the ellipsoids by changing $size
#  change the shape of the distribution by changing $offset and $spread
#  change constants in $t1 and $t2 from 1.0 to a larger number
#  to obtain elongated grains.
#######################################################################
# Z or ND 
  $test = gaussdev($i) ;
  $t3 = $size * exp(($test * $Xspread) + $Xoffset) ;

#  limit the range over which the other 2 semi-axes vary
# X or RD
  $test = gaussdev($i) ;
  $t1 = 1.0 * $t3 * exp(($test * $aspect)) ;

# Y or TD
  $test = gaussdev($i) ;
  $t2 = 4.0 * $size * exp(($test * $Yspread) + $Yoffset) ;

#######################################################################
#  Output the values to the screen
#######################################################################
  print " $t1 $t2 $t3  \n" ;
  print "1 0 0\n";
  print "0 1 0\n";
  print "0 0 1\n";
  print rand," ",rand," ",rand,"\n";
}

#
# Small ~20 grain model
#
#sub smallAxes {
#  return(0.200+rand()*0.050,0.200+rand()*0.050,0.200+rand()*0.050);
#}
#
#sub bigAxes {
#  return(0.20+rand()*0.05,0.20+rand()*0.05,0.20+rand()*0.05);
#}

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
# Large ~1000  grain model
#
sub smallAxes {
  return(0.050+rand()*0.005,0.050+rand()*0.005,0.050+rand()*0.005);
}
                                                                                
sub bigAxes {
  return(0.05+rand()*0.005,0.050+rand()*0.005,0.05+rand()*0.005);
}


sub gaussdev { 
#  remember to set $iset = 0 above

# removed the swtiching between the 2 solutions
#  that was in the NR Fortran routine
#  to avoid the setting up or passing of the flag (iset)
#  my $iset = shift(@_);
my $v1 ;
my $v2 ;
my $rsq ;
my $fac ;
#  my $gset ;
my $gasdev ;


# print STDOUT "iset = $iset \n" ;
#  if ( $iset eq 0 ) {
    do {
      $v1 = (2. * rand()) - 1.;
#      print "v1 = $v1 \n" ;
      $v2 = (2. * rand()) - 1. ;
#      print "v2 = $v2 \n" ;
      $rsq = $v1**2 + $v2**2 ;
#      print "rsq in loop = $rsq  \n" ;
    } until ( $rsq lt 1. and $rsq gt 0. ) ;
#    print  "rsq after loop = $rsq  \n" ;
    $fac = sqrt(-2. * log($rsq) / $rsq ) ;
#    $gset = $v1 * $fac ;
    $gasdev = $v2 * $fac ;
#    $iset = 1 ;
#  } {
#    $gasdev = $gset ;
#    $iset = 0 ;
#  }
# print  "result = $gasdev \n" ;
return ($gasdev) ;
}
#  ________________________________________


#  from Numerical Recipes
#      DATA iset/0/
#      if (iset.eq.0) then
#1       v1=2.*ran1(idum)-1.
#        v2=2.*ran1(idum)-1.
#        rsq=v1**2+v2**2
#        if(rsq.ge.1..or.rsq.eq.0.)goto 1
#        fac=sqrt(-2.*log(rsq)/rsq)
#        gset=v1*fac
#        gasdev=v2*fac
#        iset=1
#      else
#        gasdev=gset
#        iset=0
#      endif
#      return



