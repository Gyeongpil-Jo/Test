#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

my $col;
my $original;
my $sep = '';
my $format = "%10.5f";
my $help;
my $hb_major;

sub usage {
    exit -1;
}

GetOptions(
    "hb_major" => \$hb_major,	# include major groove hbonds
    "help" => \$help
);


my $key;
my @chains;
my $oldri;
my %curr_res;
my @curr_ch;

die unless defined $ARGV[0];

my $name = $ARGV[0];
$name =~ s/.itp//;

###################################
# Read topology
while (<>) {
    chomp;
    s/^\s+//;
    next if /^;/;

    if (/\[\s*([a-z]+)\s*\]/) {
	$key = $1;
	print STDERR "key=.$key.\n";
	next;
    }
    next unless defined $key;

    my ($ii,$at,$ri,$rn,$an,$cg,$ch,$mm) = split /\s+/;
    next unless defined $mm;

    if ($key eq "atoms") {
	if (defined $oldri) {
	    if ($ri > $oldri) {
		# add old residue to old chain
		push @curr_ch, {%curr_res};

		# start a new residue
		%curr_res = ();
		$curr_res{$an} = $ii;
		$curr_res{name} = $rn;	# resname
		$oldri = $ri;
	    }
	    elsif ($ri == $oldri) {
		# update current residue
		$curr_res{$an} = $ii;
	    }
	    else {
		# new chain
		print STDERR "starting new chain...\n";

		### flush
		if (keys %curr_res) {
		    push @curr_ch, {%curr_res};
		}

		if (@curr_ch) {
		    push @chains, [@curr_ch];
		    print STDERR "pushing new chain...\n";

		}

		# start a new residue
		%curr_res = ();
		@curr_ch = ();
		$curr_res{$an} = $ii;
		$curr_res{name} = $rn;	# resname

		$oldri = $ri;
	    }
	}
	else {
	    $oldri = $ri;
	    $curr_res{$an} = $ii;
	    $curr_res{name} = $rn;	# resname
	}

    }
    else {
    }

}

    
### flush
if (keys %curr_res) {
    push @curr_ch, {%curr_res};
}

if (@curr_ch) {
    push @chains, [@curr_ch];
    print STDERR "pushing new chain...\n";
}

###### print for test
#for my $c (0 .. $#chains) {
#    for my $r (0 .. $#{$chains[$c]}) {
#	#print "chains[$c][$r]\n";
#	for my $a (keys %{$chains[$c][$r]}) {
#	    printf "chains[%d][%d]{%s => %d}\n", $c,$r,$a,$chains[$c][$r]{$a};
#	}
#    }
#}

#### PBC bond
open FH, ">$name.pbc.itp";

print FH "; PBC bonded\n";
## PBC bond
# first P -- last O3'
print FH "[ bonds ]\n";
for my $c (0 .. $#chains) {
    my $nn = $#{$chains[$c]};
    printf FH "%8d %8d %8d\n", $chains[$c][0]{P}, $chains[$c][$nn]{"O3'"}, 1;
}

## PBC pair
# first O1P -- last C3'
# first O2P -- last C3'
# first O5' -- last C3'
# first P   -- last H3'
# first P   -- last C2'
# first P   -- last C4'
# first C5' -- last O3'
print FH "[ pairs ]\n";
for my $c (0 .. $#chains) {
    my $nn = $#{$chains[$c]};
    printf FH "%8d %8d %8d\n", $chains[$c][0]{"O1P"}, $chains[$c][$nn]{"C3'"}, 1;
    printf FH "%8d %8d %8d\n", $chains[$c][0]{"O2P"}, $chains[$c][$nn]{"C3'"}, 1;
    printf FH "%8d %8d %8d\n", $chains[$c][0]{"O5'"}, $chains[$c][$nn]{"C3'"}, 1;
    printf FH "%8d %8d %8d\n", $chains[$c][0]{"P"},   $chains[$c][$nn]{"H3'"}, 1;
    printf FH "%8d %8d %8d\n", $chains[$c][0]{"P"},   $chains[$c][$nn]{"C2'"}, 1;
    printf FH "%8d %8d %8d\n", $chains[$c][0]{"P"},   $chains[$c][$nn]{"C4'"}, 1;
    printf FH "%8d %8d %8d\n", $chains[$c][0]{"C5'"}, $chains[$c][$nn]{"O3'"}, 1;
}

## PBC angle
# last C3'    last  O3'   first P
# last O3'    first P     first O1P
# last O3'    first P     first O2P
# last O3'    first P     first O5'
print FH "[ angles ]\n";
for my $c (0 .. $#chains) {
    my $nn = $#{$chains[$c]};
    printf FH "%8d %8d %8d %8d\n",$chains[$c][$nn]{"C3'"},$chains[$c][$nn]{"O3'"},$chains[$c][0]{"P"}, 1;
    printf FH "%8d %8d %8d %8d\n",$chains[$c][$nn]{"O3'"},$chains[$c][0]{"P"},$chains[$c][0]{"O1P"}, 1;
    printf FH "%8d %8d %8d %8d\n",$chains[$c][$nn]{"O3'"},$chains[$c][0]{"P"},$chains[$c][0]{"O2P"}, 1;
    printf FH "%8d %8d %8d %8d\n",$chains[$c][$nn]{"O3'"},$chains[$c][0]{"P"},$chains[$c][0]{"O5'"}, 1;
}

## PBC pdih
# [$nn]{"C4'"}  [$nn]{"C3'"}  [$nn]{"O3'"}     [0]{"P"}  
# [$nn]{"H3'"}  [$nn]{"C3'"}  [$nn]{"O3'"}     [0]{"P"}  
# [$nn]{"C2'"}  [$nn]{"C3'"}  [$nn]{"O3'"}     [0]{"P"}  
# [$nn]{"C3'"}  [$nn]{"O3'"}      [0]{"P"}   [0]{"O1P"}
# [$nn]{"C3'"}  [$nn]{"O3'"}      [0]{"P"}   [0]{"O2P"}
# [$nn]{"C3'"}  [$nn]{"O3'"}      [0]{"P"}   [0]{"O5'"}
# [$nn]{"O3'"}      [0]{"P"}    [0]{"O5'"}   [0]{"C5'"}
print FH "[ dihedrals ]\n";
for my $c (0 .. $#chains) {
    my $nn = $#{$chains[$c]};
    printf FH "%8d %8d %8d %8d %8d\n", $chains[$c][$nn]{"C4'"}, $chains[$c][$nn]{"C3'"}, $chains[$c][$nn]{"O3'"},    $chains[$c][0]{"P"}, 9; 
    printf FH "%8d %8d %8d %8d %8d\n",$chains[$c][$nn]{"H3'"}, $chains[$c][$nn]{"C3'"}, $chains[$c][$nn]{"O3'"},    $chains[$c][0]{"P"},  9;
    printf FH "%8d %8d %8d %8d %8d\n",$chains[$c][$nn]{"C2'"}, $chains[$c][$nn]{"C3'"}, $chains[$c][$nn]{"O3'"},    $chains[$c][0]{"P"},  9;
    printf FH "%8d %8d %8d %8d %8d\n",$chains[$c][$nn]{"C3'"}, $chains[$c][$nn]{"O3'"},     $chains[$c][0]{"P"},  $chains[$c][0]{"O1P"}, 9;
    printf FH "%8d %8d %8d %8d %8d\n",$chains[$c][$nn]{"C3'"}, $chains[$c][$nn]{"O3'"},     $chains[$c][0]{"P"},  $chains[$c][0]{"O2P"}, 9;
    printf FH "%8d %8d %8d %8d %8d\n",$chains[$c][$nn]{"C3'"}, $chains[$c][$nn]{"O3'"},     $chains[$c][0]{"P"},  $chains[$c][0]{"O5'"}, 9;
    printf FH "%8d %8d %8d %8d %8d\n",$chains[$c][$nn]{"O3'"},     $chains[$c][0]{"P"},   $chains[$c][0]{"O5'"},  $chains[$c][0]{"C5'"}, 9;
}

## PBC idih

close FH;

####################################################
## hydrogen bonds restraints
open FH, ">$name.hbonds.itp";
print FH "[ bonds ]\n";
my $nn = $#{$chains[0]};
for my $i (0 .. $#{$chains[0]}) {
    my $j = $nn-$i;
    my $resi = $chains[0][$i]{name};
    my $resj = $chains[1][$j]{name};
    print STDERR "$resi($i)--$resj($j)\n";

    my $ai;
    my $aj;

    # center
    if ($resi =~ /(DA|DG|DX|RA|RG|RX)/) {	# X = 6mA
	$ai = "N1";
	$aj = "N3";
    }
    else {
	$ai = "N3";
	$aj = "N1";
    }
    printf FH "%8d %8d %8d   0.0    0.40   9.0    1000.0 ; %s(%s)--%s(%s)\n",
	    $chains[0][$i]{$ai},$chains[1][$j]{$aj},10,$resi,$ai,$resj,$aj;

    next unless defined $hb_major;

    # major
    if ($resi =~ /(DA|RA)/) {
	$ai = "N6";
	$aj = "O4";
    }
    elsif ($resi =~ /(DT|RU|DU)/) {
	$ai = "O4";
	$aj = "N6";
    }
    elsif ($resi =~ /(DG|RG)/) {
	$ai = "O6";
	$aj = "N4";
    }
    elsif ($resi =~ /(DC|RC)/) {
	$ai = "N4";
	$aj = "O6";
    }
    printf FH "%8d %8d %8d   0.0    0.50   9.0    1000.0 ; %s(%s)--%s(%s)\n",
	    $chains[0][$i]{$ai},$chains[1][$j]{$aj},10,$resi,$ai,$resj,$aj;
}
close FH;
