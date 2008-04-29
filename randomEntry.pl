######################################################################
## 
## Name: randomEntry.pl 
## Author: Steve Sintay
##
## Description: This script will read any ellipsoid file
## [ellipsoid_file] with a list of ellispoids given in the multiline
## format such as test.input or elliptical.cells and will print one
## random ellipsoid dataset to the commandline. It is intended to be
## used in place of bimodal.pl or bell_curve.pl.  
##
## usage: perl randomEntry.pl [ellipsoid_file]
##
######################################################################


# Extract the input from from the command line
$data_file = $ARGV[0];

# Open the data file and extract the data set
open(DAT, $data_file) || die("Could not open file!");
@raw_data=<DAT>;
close(DAT);

# Select a randon line from the input file
# http://www.unix.org.ua/orelly/perl/cookbook/ch08_07.htm
srand;
# SET A VARIABLE
$count = 0;

# RUN A WHILE LOOP
while ($count <= 7) {
# $rline = rand(@raw_data);

# # Determine the starting line of the dataset.
# $dline = ($rline / 6);
# $il = int ($dline);

# # Read the dataset from the file
# $line1 =  $raw_data[$il*6];
# $line2 =  $raw_data[$il*6+1];
# $line3 =  $raw_data[$il*6+2];
# $line4 =  $raw_data[$il*6+3];
# $line5 =  $raw_data[$il*6+4];
# $line6 =  $raw_data[$il*6+5];

# # Print to stdout

# # Comment out the position of the ellispsoid as we want to allow
# # recursiveSampler.pl to place it in the RVE.
# #print $line1; # Position of centroid

# print $line2; # 3 Semi axis
# print $line3; # Basis for first semi-axis
# print $line4; # Basis for second semi-axis
# print $line5; # Basis for third semi-axis
# print $line6; # Random orienation

# $count++;
my @el = randomEllipse();
my @axis = split(/\s+/,$el[0]);
my $a1 = $axis[0];
my $a = ("$a1" / 100);
my $b = ("$axis[1]" / 100);
my $c = ("$axis[2]" / 100);

print "$axis[0] $axis[1] $axis[2]\n";
print "$a $b $c\n";
print "\n";
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

# Print to stdout

# Comment out the position of the ellispsoid as we want to allow
# recursiveSampler.pl to place it in the RVE.
#print $line1; # Position of centroid

$count++;

#print $line2; # 3 Semi axis
#print $line3; # Basis for first semi-axis
#print $line4; # Basis for second semi-axis
#print $line5; # Basis for third semi-axis
#print $line6; # Random orienation

@lines = ($line2, $line3, $line4, $line5, $line6);

}
