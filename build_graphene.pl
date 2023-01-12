#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use constant PI => 4 * atan2 1, 1;

#sub numerically { $a <=> $b };

###################################
my %terminal_types = (
    0 => 'hydrogen',
    1 => 'hydroxyl',
    2 => 'carboxyl',
    3 => 'carbonyl'
);
###################################

my $help;
my $size;
my $pbc;
my $cnt;
my $hole;
GetOptions(
    "size=s" => \$size,	# e.g., "30,30" in nm
    "pbc=s"    => \$pbc, # 'xy' or 'x' or 'y'
    "cnt=s"    => \$cnt, # 'zigzag', 'zigzag-pbc', 'armchair', 'armchair-pbc'
    "hole=s"   => \$hole, # hole radius in nm
    "help"   => \$help
);

if (defined $help) {
    print "build_graphene.pl --size=XX,YY [--pbc=(xy|x|y)] [--cnt=(zigzag|armchair)]\n";
    print "\t--size=XX,YY: Create a flat graphene sheet of XX nm by YY nm.\n";
    print "\t--pbc=(xy|x|y): 'xy' puts PBC bonds to make the sheet effectively infinite. 'x' and 'y' puts PBC bonds only in one direction.\n";
    print "\t--cnt=(zigzag|zigzag-pbc|armchair|armchair-pbc): Wrap the graphene sheet to make a CNT.-pbc option makes CNT infinitely long. Chiral forms are not supported yet.\n";
    print "\t--hole=hole_radius_in_nm :.\n";
    exit
}

###################################

my $pbcx;
my $pbcy;
if (defined $pbc) {
    if ($pbc eq 'xy') {
       $pbcx = 1;
       $pbcy = 1;
    }
    elsif ($pbc eq 'x') {
       $pbcx = 1;
    }
    elsif ($pbc eq 'y') {
       $pbcy = 1;
    }
    else {
       die "'--pbc=$pbc' is illegal. See --help.\n";
    }
}

# armchair = CNT axis along x
# zigzag   = CNT axis along y
if (defined $cnt) {
    if ($cnt eq 'armchair') {
	$pbcy = 1;
	undef $pbcx;
    }
    elsif ($cnt eq 'armchair-pbc') {
	$pbcy = 1;
	$pbcx = 1;
    }
    elsif ($cnt eq 'zigzag') {
	$pbcx = 1;
	undef $pbcy;
    }
    elsif ($cnt eq 'zigzag-pbc') {
	$pbcx = 1;
	$pbcy = 1;
    }
    else {
	die "'--cnt=$cnt' is illegal. See --help.\n";
    }
}

my @dim  = split /,/, $size; # in angstrom
$dim[0] *= 10;
$dim[1] *= 10;
$dim[2] = ($dim[0] > $dim[1]) ? $dim[0] : $dim[1];
my $blen = 1.3750;  # C-C bond length in angstrom
my @ulen = (sqrt(3)*$blen, 3*$blen, 0);
my $nx = int($dim[0]/$ulen[0]);
my $ny = int($dim[1]/$ulen[1]);
my @box = ($nx*$ulen[0], $ny*$ulen[1], $dim[2]);
my @center = ($nx*$ulen[0]/2, $ny*$ulen[1]/2, $dim[2]/2);

print STDERR "A graphene layer of $dim[0] A by $dim[1] A was requested.\n";
if (defined $cnt) {
    print STDERR "A CNT will be written to graphene.pdb and graphene.itp\n";
}
else {
    print STDERR "A graphene of $box[0] by $box[1] (or $nx by $ny hexagons) will be written to graphene.pdb and graphene.itp\n";
}

if (defined $hole) {
    print STDERR "A hole of radius $hole nm will be created at the center.\n";
    $hole *= 10;    # angstrom
}
#if (defined $oterm) {
#    print STDERR "Terminal groups will be -OH.\n";
#}
#else {
#    print STDERR "Terminal groups will be -H.\n";
#}

## a building coor
#       4
#      / 
#     3
#     |
#     2
#      \
#       1
my @atom_names = ('C1','C11','O11','O12','H11',
                  'C2','C21','O21','O22','H21',
		  'C3','C31','O31','O32','H31',
		  'C4','C41','O41','O42','H41');
my @atoms;
my $rid = 0;
my $term;

# carbon (non-terminal)
for my $i (0 .. $nx-1) {
    my $x0 = $i*$ulen[0];
    for my $j (0 .. $ny-1) {
	$atoms[$i][$j]{rid} = $rid+1;

	# xy coorinates
	my $y0 = $ulen[1]*$j;
	my @xy;

	@xy = ($x0, $y0+$blen*0.5);
	if (!defined $hole || hole_dist(@xy,@center) > $hole) {
	    $atoms[$i][$j]{C2}{coor} = [ @xy ];
	    $atoms[$i][$j]{C2}{charge} = 0;
	    $atoms[$i][$j]{C2}{type} = "CA";
	}

	@xy = ($x0, $y0+$blen*1.5);
	if (!defined $hole || hole_dist(@xy,@center) > $hole) {
	    $atoms[$i][$j]{C3}{coor} = [ @xy ];
	    $atoms[$i][$j]{C3}{charge} = 0;
	    $atoms[$i][$j]{C3}{type} = "CA";
	}

	unless (!defined $pbcx && $i == ($nx-1)) {
	    # all but right terminals have C1 & C4
	    @xy = ($x0+$ulen[0]/2, $y0);
	    if (!defined $hole || hole_dist(@xy,@center) > $hole) {
		$atoms[$i][$j]{C1}{coor} = [ @xy ];
		$atoms[$i][$j]{C1}{charge} = 0;
		$atoms[$i][$j]{C1}{type} = "CA";
	    }

	    @xy = ($x0+$ulen[0]/2, $y0+$blen*2.0);
	    if (!defined $hole || hole_dist(@xy,@center) > $hole) {
		$atoms[$i][$j]{C4}{coor} = [ @xy ];
		$atoms[$i][$j]{C4}{charge} = 0;
		$atoms[$i][$j]{C4}{type} = "CA";
	    }
	}

	# bottom terminals have additional H11

	# top terminals have additional H41

	++$rid;
    }
}

#       H4    H4    H4    
#       |     |     |     
# H3    C4    C4    C4    H3   
#   \  /  \  /  \  /  \  /      
#    C3    C3    C3    C3       
#     |    |     |     |        
#    C2    C2    C2    C2       
#   /  \  /  \  /  \  /  \      
# H2    C1    C1    C1    H2   
#       |     |     |     
# H3    C4    C4    C4    H3   
#   \  /  \  /  \  /  \  /      
#    C3    C3    C3    C3       
#     |    |     |     |        
#    C2    C2    C2    C2       
#   /  \  /  \  /  \  /  \      
# H2    C1    C1    C1    H2   
#       |     |     |     
#       H1    H1    H1
##########################################
# connectivity: carbon only for now
print STDERR "Connectivity of carbon atoms.\n";
for my $i (0 .. $nx-1) {
    for my $j (0 .. $ny-1) {
	# self
	add_conn($i,$j,'C1', $i  ,$j  ,'C2');
	add_conn($i,$j,'C2', $i  ,$j  ,'C3');
	add_conn($i,$j,'C3', $i  ,$j  ,'C4');

	# right neighbor
	add_conn($i,$j,'C1', $i+1,$j  ,'C2');
	add_conn($i,$j,'C4', $i+1,$j  ,'C3');

	# bottom neighbor
	add_conn($i,$j,'C1', $i  ,$j-1,'C4');

	# left neighbor
	add_conn($i,$j,'C2', $i-1,$j  ,'C1');
	add_conn($i,$j,'C3', $i-1,$j  ,'C4');

	# top neighbor
	add_conn($i,$j,'C4', $i  ,$j+1,'C1');

    }
}

# prune a dangling carbon (carbon atoms should have 3 or 2 bonds)
for my $i (0 .. $nx-1) {
    for my $j (0 .. $ny-1) {
	for my $c ("1","2","3","4") {
	    my @cc = keys %{$atoms[$i][$j]{"C$c"}{conn}};
	    unless ($#cc == 2 || $#cc == 1) {
		print STDERR "Pruning ($i,$j,C$c) $#cc\n" if $#cc == 0;
		undef $atoms[$i][$j]{"C$c"}{coor};
	    }
	    undef %{$atoms[$i][$j]{"C$c"}{conn}};
	}
    }
}

# repeat connectivity: carbon only for now
print STDERR "Connectivity of carbon atoms.\n";
for my $i (0 .. $nx-1) {
    for my $j (0 .. $ny-1) {
	# self
	add_conn($i,$j,'C1', $i  ,$j  ,'C2');
	add_conn($i,$j,'C2', $i  ,$j  ,'C3');
	add_conn($i,$j,'C3', $i  ,$j  ,'C4');

	# right neighbor
	add_conn($i,$j,'C1', $i+1,$j  ,'C2');
	add_conn($i,$j,'C4', $i+1,$j  ,'C3');

	# bottom neighbor
	add_conn($i,$j,'C1', $i  ,$j-1,'C4');

	# left neighbor
	add_conn($i,$j,'C2', $i-1,$j  ,'C1');
	add_conn($i,$j,'C3', $i-1,$j  ,'C4');

	# top neighbor
	add_conn($i,$j,'C4', $i  ,$j+1,'C1');

    }
}

# now add terminals to carbon atoms with connectivity less than 3.
for my $i (0 .. $nx-1) {
    for my $j (0 .. $ny-1) {
	#for my $k (keys %{$atoms[$i][$j]{C1}{conn}}) {
	#    print STDERR "$atoms[$i][$j]{C1}{conn}{$k}{coor}[0]\n";
	#}
	
	# terminals
	for my $c ("1","2","3","4") {
	    next unless defined $atoms[$i][$j]{"C$c"}{coor};
	    my @cc = keys %{$atoms[$i][$j]{"C$c"}{conn}};
	    next if $#cc == 2;

	    my $term = determine_terminal();
            $atoms[$i][$j]{"C$c"}{term} = $term;
	    #print STDERR "Add $term terminal to ($i,$j,C$c) = $#cc\n";

	    my @ccc = @{$atoms[$i][$j]{"C$c"}{coor}};
	    my @cc0;
	    my $cc0n = 0;
	    # neighbor carbons connnected to this carbon
	    for my $k (keys %{$atoms[$i][$j]{"C$c"}{conn}}) {
		$cc0[0] += $atoms[$i][$j]{"C$c"}{conn}{$k}{coor}[0];
		$cc0[1] += $atoms[$i][$j]{"C$c"}{conn}{$k}{coor}[1];
		++$cc0n;
	    }
	    #print STDERR "cc0n = $cc0n\n";
	    my @dx = ($ccc[0] - $cc0[0]/$cc0n, $ccc[1] - $cc0[1]/$cc0n);
	    $dx[0] *= 3;
	    $dx[1] *= 3;

	    if ($term eq 'hydrogen') {
		$atoms[$i][$j]{"C$c"}{charge} = -0.115;

		# add hydrogen
		$atoms[$i][$j]{"H$c"."1"}{coor} = [ $ccc[0]+$dx[0], $ccc[1]+$dx[1]];
                $atoms[$i][$j]{"H$c"."1"}{charge} = 0.115;
                $atoms[$i][$j]{"H$c"."1"}{type} = "HA";
		add_conn($i,$j,"C$c", $i  ,$j,"H$c"."1");
	    }
	    elsif ($term eq 'hydroxyl') {
		$atoms[$i][$j]{"C$c"}{charge} = 0.25;

		# add oxygen
		$atoms[$i][$j]{"O$c"."1"}{coor} = [ $ccc[0]+$dx[0], $ccc[1]+$dx[1]];
                $atoms[$i][$j]{"O$c"."1"}{charge} = -0.65;
                $atoms[$i][$j]{"O$c"."1"}{type} = "OH";
		add_conn($i,$j,"C$c", $i  ,$j,"O$c"."1");

		# add hydrogen
		$atoms[$i][$j]{"H$c"."1"}{coor} = [ $ccc[0]+$dx[0]*2, $ccc[1]+$dx[1]*2];
                $atoms[$i][$j]{"H$c"."1"}{charge} = 0.4;
                $atoms[$i][$j]{"H$c"."1"}{type} = "HO";
		add_conn($i,$j,"O$c"."1", $i  ,$j,"H$c"."1");
	    }
	    elsif ($term eq 'carboxyl') {
		$atoms[$i][$j]{"C$c"}{charge} = 0.0;

		# add carbon
		$atoms[$i][$j]{"C$c"."1"}{coor} = [ $ccc[0]+$dx[0], $ccc[1]+$dx[1]];
                $atoms[$i][$j]{"C$c"."1"}{charge} = 0.8;
                $atoms[$i][$j]{"C$c"."1"}{type} = "C";
		add_conn($i,$j,"C$c", $i  ,$j,"C$c"."1");

		# add oxygen1
		$atoms[$i][$j]{"O$c"."1"}{coor} = [ $ccc[0]+$dx[0]*2, $ccc[1]+$dx[1]*2];
                $atoms[$i][$j]{"O$c"."1"}{charge} = -0.9;
                $atoms[$i][$j]{"O$c"."1"}{type} = "O2";
		add_conn($i,$j,"C$c"."1", $i  ,$j,"O$c"."1");

		# add oxygen2
		$atoms[$i][$j]{"O$c"."2"}{coor} = [ $ccc[0]+$dx[0]*2, $ccc[1]+$dx[1]*2.5];
                $atoms[$i][$j]{"O$c"."2"}{charge} = -0.9;
                $atoms[$i][$j]{"O$c"."2"}{type} = "O2";
		add_conn($i,$j,"C$c"."1", $i  ,$j,"O$c"."2");
	    }
	    elsif ($term eq 'carbonyl') {
	    }
	}
    }
}

#######################################
# assign atom id
my $aid = 0;
for my $i (0 .. $nx-1) {
    for my $j (0 .. $ny-1) {
	for my $a (@atom_names) {
	    if (defined $atoms[$i][$j]{$a}{coor}) {
		$atoms[$i][$j]{$a}{aid} = ++$aid;
	    }
	}
    }
}
my $natoms = $aid;


# build connectivity map
my @conn;
for my $i (0 .. $nx-1) {
    for my $j (0 .. $ny-1) {
	for my $a (@atom_names) {
	    if (defined $atoms[$i][$j]{$a}{conn}) {
		my $aid1 = $atoms[$i][$j]{$a}{aid};
		for my $k (keys %{$atoms[$i][$j]{$a}{conn}}) {
		    #    print STDERR "$atoms[$i][$j]{C1}{conn}{$k}{coor}[0]\n";
		    my $aid2 = $atoms[$i][$j]{$a}{conn}{$k}{aid};
		    push @{$conn[$aid1]}, $aid2;
		    push @{$conn[$aid2]}, $aid1;
		}
	    }
	}
    }
}

###############################################################
# print PDB
print STDERR "Print PDB...\n";
open FPDB, ">graphene.pdb";
#CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1
printf FPDB "CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00 P 1           1\n",
	@box;
#ATOM      1  OW  SOL     1       5.690  12.751  11.651  1.00  0.00
my $format = "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n";
for my $i (0 .. $nx-1) {
    for my $j (0 .. $ny-1) {
	for my $a (@atom_names) {
	    next unless defined $atoms[$i][$j]{$a}{aid};
	    my @xyz = (@{$atoms[$i][$j]{$a}{coor}}, 0.5*$box[2]);
	    if (defined $cnt && ( $cnt eq 'zigzag' or $cnt eq 'zigzag-pbc')) {
		# CNT axis = y, wrap about x.
		my $R = $box[0]/(2*PI);
		my $theta = $xyz[0]/$R;
		$xyz[0] = $R*sin($theta);
		$xyz[2] = $R*(1-cos($theta));
	    }
	    if (defined $cnt && ( $cnt eq 'armchair' or $cnt eq 'armchair-pbc')) {
		# CNT axis = x, wrap about y.
		my $R = $box[1]/(2*PI);
		my $theta = $xyz[1]/$R;
		$xyz[1] = $R*sin($theta);
		$xyz[2] = $R*(1-cos($theta));
	    }
	    my @ccc = ( "ATOM", 
			$atoms[$i][$j]{$a}{aid} % 100000, 
			$a, '', "GRA", '',
			$atoms[$i][$j]{rid} % 10000, 
			'',
			@xyz, 0.0, 0.0);
	    printf FPDB $format, @ccc;
	}
    }
}

################################################################
## print top
print STDERR "Print Gromacs topology...\n";
open FTOP, ">graphene.itp";
open FRES, ">graphene.posres.itp";
print FTOP ";build_graphene.pl nx=$nx ny=$ny\n";
print FTOP "[ moleculetype ]\n";
print FRES "[ position_restraints ]\n";
print FRES "; ai  funct  fcx    fcy    fcz\n";
print FTOP "; Name            nrexcl\n";
print FTOP "graphene           3\n\n";

print FTOP "[ atoms ]\n";
print FTOP ";  nr  type  resnr residue atom   cgnr charge  mass\n";

for my $i (0 .. $nx-1) {
    for my $j (0 .. $ny-1) {
	for my $a (@atom_names) {
	    next unless defined $atoms[$i][$j]{$a}{aid};
	    my $charge = $atoms[$i][$j]{$a}{charge};
	    my $type = $atoms[$i][$j]{$a}{type};
	    my $atom;
	    my $mass;
	    if ($a =~ /C/) {
		$mass = 12.011;
		$atom = "C$i";
	    }
	    if ($a =~ /O/) {
		$mass = 16.000;
		$atom = "O$i";
	    }
	    if ($a =~ /H/) {
		$mass = 1.008;
		$atom = "H$i";
	    }

	    printf FTOP "%4d%7s%8d%7s%7s%8d%8.3f%8.3f\n",
		    $atoms[$i][$j]{$a}{aid},    # nr
		    $type,		 # type
		    $atoms[$i][$j]{rid},        # resnr
		    #"GRA",                      # residue
		    #$a,                         # atom
		    ##### use residue and atom names to indicate colorum, row.
		    "G$j",			# residue
		    $atom,		      # atom
		    $atoms[$i][$j]{$a}{aid},    #$i,            # cgnr
		    $charge,                    # charge
		    $mass;			# mass

	    printf FRES "%4d%8d%8d%8d%8d\n", $atoms[$i][$j]{$a}{aid}, 1, 1000, 1000, 1000 if ($type eq "CA");
	}
    }
}

## build angles, dihedrals;
my %bonds;
my %angles;
my %dihedrals;
my %pairs;
for my $a1 (1 .. $#conn) {
    next unless defined $conn[$a1];
    printf FPDB "CONECT%5d", $a1;

    for my $a2 (@{$conn[$a1]}) {
	my @tmp = ($a1 > $a2) ? ($a2,$a1,1) : ($a1,$a2,1);
	$bonds{sprintf "%8d%8d%8d\n", @tmp} = 1;

	# PDB CONECT
	printf FPDB "%5d", $a2;

	next unless defined $conn[$a2];
	for my $a3 (@{$conn[$a2]}) {
	    next if $a1 == $a3;
	    #my @tmp = ($a1 > $a3) ? ($a3,$a2,$a1,5) : ($a1,$a2,$a3,5); urey-bradley
	    my @tmp = ($a1 > $a3) ? ($a3,$a2,$a1,1) : ($a1,$a2,$a3,1);
	    $angles{sprintf "%8d%8d%8d%8d\n", @tmp} = 1;
	    next unless defined $conn[$a3];
	    for my $a4 (@{$conn[$a3]}) {
		next if $a1 == $a4;
		next if $a2 == $a4;
		my @tmp = ($a1 > $a4)?($a4,$a3,$a2,$a1,9):($a1,$a2,$a3,$a4,9);
		$dihedrals{sprintf "%8d%8d%8d%8d%8d\n", @tmp} = 1;
		$pairs{sprintf "%8d%8d%8d\n", $tmp[0],$tmp[3],1} = 1;
	    }
	}
    }
    printf FPDB "\n";
}

print FTOP "\n[ bonds ]\n";
print FTOP ";  ai    aj funct            c0            C1            C2            C3\n";
print FTOP sort keys %bonds;

print FTOP "\n[ pairs ]\n";
print FTOP ";  ai    aj funct            c0            C1            C2            C3\n";
print FTOP sort keys %pairs;

print FTOP "\n[ angles ]\n";
print FTOP ";  i    aj    ak funct            c0            C1            C2            C3\n";
print FTOP sort keys %angles;

print FTOP "\n[ dihedrals ]\n";
print FTOP ";  ai    aj    ak    al funct            c0            C1            C2            C3            C4            c5\n";
print FTOP sort keys %dihedrals;

close FTOP;
close FRES;
close FPDB;

exit 0;

sub add_conn {
    my $i1 = shift;
    my $j1 = shift;
    my $a1 = shift;
    my $i2 = shift;
    my $j2 = shift;
    my $a2 = shift;

    # pbc
    $i1 -= $nx if ($i1 >= $nx && defined $pbcx);
    $i1 += $nx if ($i1 < 0    && defined $pbcx);
    $j1 -= $ny if ($j1 >= $ny && defined $pbcy);
    $j1 += $ny if ($j1 < 0    && defined $pbcy);
    $i2 -= $nx if ($i2 >= $nx && defined $pbcx);
    $i2 += $nx if ($i2 < 0    && defined $pbcx);
    $j2 -= $ny if ($j2 >= $ny && defined $pbcy);
    $j2 += $ny if ($j2 < 0    && defined $pbcy);

    return if $i1 < 0 || $i1 >= $nx;
    return if $j1 < 0 || $j1 >= $ny;
    return if $i2 < 0 || $i2 >= $nx;
    return if $j2 < 0 || $j2 >= $ny;

    if (defined $atoms[$i1][$j1]{$a1}{coor} && defined $atoms[$i2][$j2]{$a2}{coor}) {
	$atoms[$i1][$j1]{$a1}{conn}{"$i2-$j2-$a2"} = $atoms[$i2][$j2]{$a2};
	$atoms[$i2][$j2]{$a2}{conn}{"$i1-$j1-$a1"} = $atoms[$i1][$j1]{$a1};
    }
#    if (defined $atoms[$i1][$j1]{$a1}{aid} && defined $atoms[$i2][$j2]{$a2}{aid}) {
#	#$conn_mat[$atoms[$i1][$j1]{$a1}{aid}][$atoms[$i2][$j2]{$a2}{aid}] = 1;
#	#$conn_mat[$atoms[$i2][$j2]{$a2}{aid}][$atoms[$i1][$j1]{$a1}{aid}] = 1;
#	push @{$conn[$atoms[$i1][$j1]{$a1}{aid}]}, $atoms[$i2][$j2]{$a2}{aid};
#	push @{$conn[$atoms[$i2][$j2]{$a2}{aid}]}, $atoms[$i1][$j1]{$a1}{aid};
#    }
#    else {
#	print STDERR "Trying to connect undefined atoms: ($i1,$j1,$a1) - ($i2,$j2,$a2)\n";
#    }
}

sub determine_terminal {
    my $t = int(rand(3));
    #print STDERR "Use $terminal_types{$t}\n";
    return $terminal_types{$t};
}

sub hole_dist {
    my @r = ($_[0]-$_[2],$_[1]-$_[3]);
    return sqrt($r[0]*$r[0]+$r[1]*$r[1]);
}
