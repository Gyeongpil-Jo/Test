#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

###################################
#my $orig_pattern = "(ATOM  )(.{5}).(.{4})(.)(.{3}).(.)(.{4})(.)...(.{8})(.{8})(.{8})(.{6})(.{6})";
my $chm_pattern = "(ATOM  )(.{5}).(.{4})(.)(.{4})(.)(.{4})(.)...(.{8})(.{8})(.{8})(.{6})(.{6})";
#
#my $orig_format = "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n"
my $chm_format = "%-6s%5d %-4s%1s%4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n";
#
#my @orig_len = (6,5,4,1,3,1,4,1,8,8,8,6,6);
my @chm_len  = (6,5,4,1,4,1,4,1,8,8,8,6,6);
#
####################################
## Default
my $pattern = $chm_pattern;
my $format = $chm_format;
my @len = @chm_len;
####################################
# 
###################################
my $help;
my $phi = 0;
GetOptions(
    "phi=s" => \$phi,	# degree per angstrom; + is untwisting. -0.7 for PBC DNA for NAB output?
    "help" => \$help
);

$phi = $phi * 3.14 / 180;
###################################

if (defined($help)) {
    print STDERR "\n";

    exit;
}

my @data;
my @CM;
my $nn;
while (<>) {
    print if /CRYST/;
    chomp;

    next unless /$pattern/;

    my ($atom,$x,$y,$z) = ($3,$9,$10,$11);
    my @tmp = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13);
    push @data, [@tmp];

    # compute CM
    unless ($atom =~ /H/) {
	$CM[0] += $x;
	$CM[1] += $y;
	$CM[2] += $z;
	++$nn;
    }
}
$CM[0] /= $nn;
$CM[1] /= $nn;
$CM[2] /= $nn;
print STDERR "CM = (@CM)\n";

for my $p (@data) {
    my ($x,$y,$z) = ($p->[8]-$CM[0],$p->[9]-$CM[1],$p->[10]-$CM[2]);
    my $theta = $phi*$z;
    $p->[8]  = $CM[0] + $x * cos($theta) - $y * sin($theta);
    $p->[9]  = $CM[1] + $x * sin($theta) + $y * cos($theta);
    $p->[10] = $CM[2] + $z;

    printf $format, @{$p};
}
print "END\n";

exit;
