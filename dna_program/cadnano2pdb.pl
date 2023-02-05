#!/usr/bin/perl -w

###################################
# We need JSON.pm http://search.cpan.org/~makamaka/JSON-2.53/lib/JSON.pm
# Install JSON.pm:
# - perl Makefile.PL PREFIX=$HOME/myusr/perl5
# - make
# - make test
# - make install
# Then, put the following line in .bashrc or .profile_user etc.
# - export PERL5LIB=$PERL5LIB:$HOME/myusr/perl5/lib64/perl5:$HOME/myusr/perl5/share/perl5
# You might need to install ExtUtils-MakeMaker http://search.cpan.org/~mschwern/ExtUtils-MakeMaker-6.64/lib/ExtUtils/MakeMaker.pm
#
###################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
#use JSON::PP;
use JSON;
use Cwd 'abs_path';
my $PLPATH = dirname(abs_path($0));

###################################
my $lattice;
#my $chainsep;
my $noloop;
my $noskip;
my $forcefield;
my $scafseqname;
my $help;
my $csvfile;
my $scaf_resi=0;    # global variable for counting multiple scaffold strands.
GetOptions(
    "ff=s" => \$forcefield,
    "lattice=s" => \$lattice,
    "scaf=s" => \$scafseqname,
    "stap=s" => \$csvfile,
    #"chainsep" => \$chainsep,
    "ignore-loop" => \$noloop,
    "ignore-skip" => \$noskip,
    "help" => \$help
);

###################################
my $dd = 23.0;	# interaxial distance in A

#print ord('A');	=> 65
#print ord('Z');	=> 90

####################################
# Read sequence files
my %scafseq;

my @seq_files = `ls $PLPATH/cadnano2pdb.seq.*.dat`;

for my $f (@seq_files) {
    open FH, $f;
    $f =~ /.*cadnano2pdb.seq.(.*).dat/;
    my $name = $1;
    #print STDERR "Reading sequence file: \"$name\"\n";
    my @tmp;
    while (<FH>) {
	chomp;
	if (/^;(.*)/) {
	    #print STDERR "\"$name\": $1\n";
	    next;
	}
	s/[^atgcATGC]//g;
	tr/atgc/ATGC/;
	push @tmp, split //;
    }
    $scafseq{$name} = [@tmp];
    close FH;
}
#$scafseq{'none'} = [];
unless (%scafseq) {
    die "Error: no sequence files are found. Place sequence files in the same directory as the perl script.\n";
}

if (defined($help) or ($#ARGV < 0)) {
    &print_help;
    exit;
}

my $basename = [split /\./,basename($ARGV[0])]->[0];
print STDERR "json=$basename.json\n";
#print STDERR "stap=$basename.csv\n";

### select lattice type.
if (defined $lattice && !( $lattice eq "honeycomb" or $lattice eq "square") ) {
    print STDERR "Warning: --lattice=$lattice is undefined. Choose again!\n";
    undef $lattice;
}

unless (defined $lattice) {
    do {
	print STDERR "Select lattice type: honeycomb or square.\n";
	chomp($lattice = <STDIN>);
    }
    until ($lattice eq "honeycomb" or $lattice eq "square");
}
print STDERR "Lattice type \"$lattice\" is selected.\n";

### select forcefield type.
if (defined $forcefield && !( $forcefield eq "charmm" or $forcefield eq "amber") ) {
    print STDERR "Warning: --forcefield=$forcefield is undefined. Choose again!\n";
    undef $forcefield;
}
unless (defined $forcefield) {
    do {
	print STDERR "Select forcefield type: charmm or amber.\n";
	chomp($forcefield = <STDIN>);
    }
    until ($forcefield eq "charmm" or $forcefield eq "amber");
}
print STDERR "Forcefield type \"$forcefield\" is selected.\n";

#unless (defined $lattice && ($lattice eq "honeycomb" or $lattice eq "square") 
#	&& defined $forcefield && ($forcefield eq "charmm" or $forcefield eq "amber")) {
#    print STDERR "json2pdb.pl --lattice=(honeycomb|square) --ff=(charmm|amber)\n";
#    exit;
#}

####################################
# Read PDB of A,T,G,C
# - AT and GC pairs by rotating one NT by 180 wrt x axis.
# - AT and GC pairs are parallel to x axis.
my %pdb;

if ($forcefield eq "charmm") {
    print STDERR "CHARMM36 force field\n";
    &init_pdb_charmm;
}
else {
    print STDERR "AMBER99 force field\n";
    &init_pdb_amber;
}

# determine sequence
if (defined $scafseqname && !defined($scafseq{$scafseqname}) ) {
    print STDERR "Warning: --scaf=$scafseqname is not defined.\n";
    print STDERR "\t\tSCAFFOLD_SEQUENCE=\n";
    for my $s (keys %scafseq) {
	printf STDERR "\t\t\t %-15s (%5d nts)\n", $s, $#{$scafseq{$s}}+1;
    }
    #undef $scafseqname;
    exit;
}

if (!defined $scafseqname && !defined $csvfile) {
    print STDERR "Warning: at least one of --scaf and --stap options must be given\n";
    exit;
}

#unless (defined $scafseqname) {
#    do {
#    print STDERR "----------------------------\n";
#    print STDERR "Select scaffold sequence among these:\n";
#	for my $s (keys %scafseq) {
#	    printf STDERR ">> %-15s (%5d nts)\n", $s, $#{$seq{$s}}+1;
#	}
#	chomp($scafseqname = <STDIN>);
#    }
#    until (defined($scafseq{$scafseqname}));
#    print STDERR "\"$scafseqname\" is selected\n";
#    print STDERR "----------------------------\n";
#}

####################################
# Ideal ENM distances
my %enm_stack;
my %enm_bp;
&init_enm;
my @baseatoms = (' N1 ',' C2 ',' N3 ',' C4 ',' C5 ',' C6 ',' N7 ',' C8 ',' N9 ');

####################################
# Read sequence (CSV)
my @csv;
if (defined $csvfile) {
    open FH, $csvfile or die "Can't open $csvfile\n";
    print STDERR "Reading CVS file $csvfile\n";
    while (<FH>) {
	next unless /([0-9]+)\[([0-9]+)\],([0-9]+)\[([0-9]+)\],([ATGC]+),[0-9]+/;

	my %tmp = (
	    '5helix' => $1,
	    '5base' => $2,
	    '3helix' => $3,
	    '3base' => $4
	);
	$tmp{seq} = [split('',$5)];
	push @csv, {%tmp};	# {} for hash reference, [] for array reference.
	#print STDERR "$1 $2 $3 $4 $5\n";
    }
    close FH;

    print STDERR "#staples in $csvfile = $#csv+1\n";
    #print STDERR $csv[0]{seq}[0]."\n";
}


####################################
# Read json
# Indices for helix and base starts from 0.
# @vstrands array has a hash of these:
#	num : helix index from 0. NOTE that @vstrands order and $num don't match! 
#	      Use ih2id() function.
#	      Order of helix id doesn't matter because helix positions determined by (row,col)
#	stap	: connectivity array of bases [5' helix, 5' base, 3' helix, 3' base]
#	scaf	: connectivity array of bases
#	row	
#	col	
#	stap_colors	
#	skip	    ; deletion
#	loop	    ; insertion
#	scafLoop    ; deprecated?
#	stapLoop    ; deprecated?
### defined by me #############################
#	scaf_nt    : sequence
#	stap_nt    : sequence
#	stap_ri    : staple residue id (one-based)
#	scaf_ri    : scaffold residue id (one-based)
#	theta	    : rotation considering insertion/deletion.
#	z	    : z coordinate considering insertion/deletion.
#	stap_atoms : hash of atom name & atom index
#	scaf_atoms : hash of atom name & atom index
#	hb_major	: two atom indices of H-bond in major groove. DON'T USE IT!
#	hb_minor	: two atom indices of H-bond in minor groove. DON'T USE IT!
#	hb_center	: two atom indices of H-bond in center. DON'T USE IT!
#	c1p_scaf		: atom index of C1' atom DON'T USE IT!
#	c1p_stap		: atom index of C1' atom DON'T USE IT!
#	( N1 | C2 | N3 | C4 | C5 | C6 | N7 | C8 | N9 ) : base atom indices
####################################

my $json = decode_json(<>);
my @vstrands = @{$json->{"vstrands"}};
## build ih2id array.
my @ih2id;
for my $i (0 .. $#vstrands) {
    my $ih = $vstrands[$i]->{num};
    $ih2id[$ih] = $i;
}

####################### 
# renumber $ih so that we can forget about it.
for my $ih (0 .. $#vstrands) {
    # scaffold
    for my $ib (0 .. $#{$vstrands[$ih]->{scaf}}) {
	my ($pih,$pib,$nih,$nib) = @{$vstrands[$ih]{scaf}[$ib]};
	$vstrands[$ih]{scaf}[$ib][0] = $ih2id[$pih] if ($pih > -1);
	$vstrands[$ih]{scaf}[$ib][2] = $ih2id[$nih] if ($nih > -1);
    }
    # staples
    for my $ib (0 .. $#{$vstrands[$ih]->{stap}}) {
	my ($pih,$pib,$nih,$nib) = @{$vstrands[$ih]{stap}[$ib]};
	$vstrands[$ih]{stap}[$ib][0] = $ih2id[$pih] if ($pih > -1);
	$vstrands[$ih]{stap}[$ib][2] = $ih2id[$nih] if ($nih > -1);
    }
}
for my $i (0 .. $#csv) {
    $csv[$i]{'5helix'} = $ih2id[$csv[$i]{'5helix'}];
}

my $global_ri = 0;
#for my $ih (0 .. $#vstrands) {
#    print STDERR "stap[$ih]: $#{$vstrands[$ih]->{stap}}\n";
#    print STDERR "scaf[$ih]: $#{$vstrands[$ih]->{scaf}}\n";
#    print STDERR "loop[$ih]: $#{$vstrands[$ih]->{loop}}\n";
#    print STDERR "skip[$ih]: $#{$vstrands[$ih]->{skip}}\n";
#}
### Ignore skip or loop if needed. csv should be also modified accordingly.

if (defined $noloop) {
    for my $ih (0 .. $#vstrands) {
	for my $ib (0 .. @{$vstrands[$ih]->{loop}}) {
	    $vstrands[$ih]->{loop}[$ib] = 0;
	}
    }
}
if (defined $noskip) {
    for my $ih (0 .. $#vstrands) {
	for my $ib (0 .. @{$vstrands[$ih]->{skip}}) {
	    $vstrands[$ih]->{skip}[$ib] = 0;
	}
    }
}

### example
# @{$vstrands[0]{stap}[10]} has four integers of [5' helix, 5' base, 3' helix, 3' base] for helix 0 and base 10.
#printf STDERR "%d %d %d %d\n", @{$vstrands[0]{stap}[10]};

## initialize sequence array
#for my $vs (@vstrands) {
#    $vs->{stap_nt} = [];
#}

####################################################
# determine theta and z considering insertion/deletion
# For honeycomb, 2 turns (720 degrees) per 21 bps.
# For square, 3 turns (1080 degrees) per 32 bps.
my $nblock = ($lattice eq "honeycomb") ? 7 : 8;
my $theta_per_bp = ($lattice eq "honeycomb") ? 360*2/21 : 360*3/32;
my $z_per_bp = ($lattice eq "honeycomb") ? 3.4 : 3.4;
for my $ih (0 .. $#vstrands) {
    my $ib = 0;
    do {
	my $dn = 0;	# change
	for my $i ($ib .. ($ib+$nblock-1)) {
	    $dn += $vstrands[$ih]->{loop}[$i]+$vstrands[$ih]->{skip}[$i];
	}

	my $resi=0;
	my $theta0 = (int($ib/$nblock))*$nblock*$theta_per_bp;
	my $z0 = (int($ib/$nblock))*$nblock*$z_per_bp;
	for my $i ($ib .. ($ib+$nblock-1)) {
	    next if $vstrands[$ih]->{skip}[$i] == -1;

	    #########################################
	    # important convention for loops:
	    # For upstream (+z direction), theta & z go like 0, 1, 2 ...
	    # For downstream (-z direction), theta & z go in opposite way.
	    for my $il (0 .. $vstrands[$ih]->{loop}[$i]) {
		my $theta = $theta0 +$theta_per_bp*$resi*$nblock/($nblock+$dn);
		my $z = $z0+$z_per_bp*$resi*$nblock/($nblock+$dn);
		push @{$vstrands[$ih]->{theta}[$i]},$theta;
		push @{$vstrands[$ih]->{z}[$i]},$z;

		++$resi;
		#print STDERR "($ih,$i) dn=$dn theta=$theta z=$z\n";

		# initialize H-bond indices
		@{$vstrands[$ih]{hb_major}[$i][$il]} = ();
		@{$vstrands[$ih]{hb_center}[$i][$il]} = ();
		@{$vstrands[$ih]{hb_minor}[$i][$il]} = ();
		#$vstrands[$ih]{c1p_stap}[$i][$il];
		#$vstrands[$ih]{c1p_scaf}[$i][$il];
	    }
	}

	$ib += $nblock;
    }
    while ($ib < $#{$vstrands[$ih]->{stap}});
}

my $ch = 65;	# chain starting from A
my $nseg = 1;

####################################################
# Assign sequences to scaffolds & staples 
print STDERR "Start assigning sequences............\n";
for my $ih (0 .. $#vstrands) {
    for my $ib (0 .. $#{$vstrands[$ih]->{stap}}) {
	if (defined $csvfile) {
	    my ($pih,$pib,$nih,$nib) = @{$vstrands[$ih]{stap}[$ib]};
	    if ($pib == -1 && $nib != -1) { # found 5' end
		assign_seq_staple($ih,$ib);
	    }
	}
	if (defined $scafseqname) {
	    my ($pih,$pib,$nih,$nib) = @{$vstrands[$ih]{scaf}[$ib]};
	    if ($pib == -1 && $nib != -1) { # found 5' end
		assign_seq_scaffold($ih,$ib);
	    }
	}
    }
}

####################################################
# Check sequence complementarity
print STDERR "Checking sequence complementarity............\n";
#printf STDERR "ib=    ";
#for my $ib (0 .. $#{$vstrands[0]->{stap}}) {
#    print STDERR "|" if (($ib%8)==0);
#    print STDERR $ib%$nblock;
#}
#print STDERR "\n";
for my $ih (0 .. $#vstrands) {
    #printf STDERR "ih=%3d ", $ih;
    for my $ib (0 .. $#{$vstrands[$ih]->{stap}}) {
	#print STDERR "|" if (($ib%8)==0);
	my @list_il = (0 .. $vstrands[$ih]{loop}[$ib]);
	for my $il (@list_il) {
	    if (defined $vstrands[$ih]{stap_nt}[$ib][$il]) {
		if (defined $vstrands[$ih]{scaf_nt}[$ib][$il]) {
		    print STDERR "Warning: scaf=$vstrands[$ih]{scaf_nt}[$ib][$il] stap=$vstrands[$ih]{stap_nt}[$ib][$il] at ih=$ih and ib=$ib\n" 
		    unless ($vstrands[$ih]{stap_nt}[$ib][$il] eq wc($vstrands[$ih]{scaf_nt}[$ib][$il]));
		}
		else {
		    $vstrands[$ih]{scaf_nt}[$ib][$il] = wc($vstrands[$ih]{stap_nt}[$ib][$il]);
		}
	    }
	    elsif (defined $vstrands[$ih]{scaf_nt}[$ib][$il]) {
		$vstrands[$ih]{stap_nt}[$ib][$il] = wc($vstrands[$ih]{scaf_nt}[$ib][$il]);
	    }
	    #print STDERR $scaf_nt;
	}
    }
    #print STDERR "\n";
}

####################################################
# Assign sequences to scaffolds & staples 
print STDERR "Checking scaffold sequence............\n";
for my $ih (0 .. $#vstrands) {
    for my $ib (0 .. $#{$vstrands[$ih]->{scaf}}) {
	my ($pih,$pib,$nih,$nib) = @{$vstrands[$ih]{scaf}[$ib]};
	if ($pib == -1 && $nib != -1) { # found 5' end
	    print STDERR "Found a 5'-end of scaffold at ih=$ih,ib=$ib\n";
	    check_seq_scaf($ih,$ib);
	}
    }
}
print STDERR "Checking staple sequence............\n";
for my $ih (0 .. $#vstrands) {
    for my $ib (0 .. $#{$vstrands[$ih]->{stap}}) {
	my ($pih,$pib,$nih,$nib) = @{$vstrands[$ih]{stap}[$ib]};
	if ($pib == -1 && $nib != -1) { # found 5' end
	    print STDERR "Found a 5'-end of staple at ih=$ih,ib=$ib\n";
	    check_seq_stap($ih,$ib);
	}
    }
}

#################################################
# Print scaffold sequence for checking
printf STDERR "ib=    ";
for my $ib (0 .. $#{$vstrands[0]->{scaf}}) {
    print STDERR "|" if (($ib%8)==0);
    print STDERR $ib%$nblock;
}
print STDERR "\n";
for my $ih (0 .. $#vstrands) {
    printf STDERR "ih=%3d ", $ih;
    for my $ib (0 .. $#{$vstrands[$ih]->{scaf}}) {
	print STDERR "|" if (($ib%8)==0);
	my @list_il = (0 .. $vstrands[$ih]{loop}[$ib]);
	for my $il (@list_il) {
	    print STDERR (defined $vstrands[$ih]{scaf_nt}[$ib][$il]) ?
		$vstrands[$ih]{scaf_nt}[$ib][$il] : ".";
	}
    }
    print STDERR "\n";
}

####################################################
# Print scaffold
# now supports multiple scaffolds.
print STDERR "Printing scaffold PDB............\n";
for my $ih (0 .. $#vstrands) {
    for my $ib (0 .. $#{$vstrands[$ih]->{scaf}}) {
	my ($pih,$pib,$nih,$nib) = @{$vstrands[$ih]{scaf}[$ib]};
	if ($pib == -1 && $nib != -1) { # found 5' end
	    #print STDERR "Found a scaffold at [$ih,$ib] = [$pih,$pib,$nih,$nib]\n";
	    trace_scaffold($ih,$ib);
	    ++$nseg;
	}
    }
}

####################################################
# Trace staples 
my $global_aid=0;	# atom index
print STDERR "Printing staple PDB............\n";
for my $ih (0 .. $#vstrands) {
    for my $ib (0 .. $#{$vstrands[$ih]->{stap}}) {
	#printf STDERR "@{$base}\n";

	my ($pih,$pib,$nih,$nib) = @{$vstrands[$ih]{stap}[$ib]};
	if ($pib == -1 && $nib != -1) { # found 5' end
	    #print STDERR "Found a staple at [$ih,$ib] = [$pih,$pib,$nih,$nib]\n";
	    trace_staple($ih,$ib);
	    ++$nseg;
	}
    }
}

####################################################
# inter-helical push 
# usage:
# > cat push.helix.for.make_ndx  | make_ndx  -f hextube.pdb -o tmp.ndx
# > awk '/^[0-9]/ && NF == 2 {printf "bond %8d %8d %8d %8d\n", $1-1, $2-1, 1, 30}' tmp.ndx > push.extrabonds
open FH, ">push.helix.for.make_ndx";
print FH "del 0 - 30\n";
my @nnih; # nearest neighbor list of helices
for my $ih (0 .. $#vstrands) {
    my $dir = is_upstream($ih); # +1 for upstream, -1 for downstream.
    for my $ib (0 .. $#{$vstrands[$ih]->{stap}}) {
	next if ($vstrands[$ih]{stap}[$ib][1] == -1 
	      && $vstrands[$ih]{stap}[$ib][3] == -1);
	next if $vstrands[$ih]->{skip}[$ib] == -1;

	my ($stap_pih,$stap_pib,$stap_nih,$stap_nib) = @{$vstrands[$ih]{stap}[$ib]};
	my ($scaf_pih,$scaf_pib,$scaf_nih,$scaf_nib) = @{$vstrands[$ih]{scaf}[$ib]};

	if ($ih >= 0 && $stap_nih >= 0 && $ih != $stap_nih) {
	    #print STDERR "stap cross from H$ih to H$stap_nih at $ib\n";
	    # NOTE. I'm not sure how many base pairs from HJ need to be excluded.
	    # Obviously, must be lattice dependent.
	    for my $off (0 .. 7) {
		$nnih[$ih][$stap_nih][$ib+$off] = 1;
		$nnih[$stap_nih][$ih][$ib+$off] = 1;
		$nnih[$ih][$stap_nih][$ib-$off] = 1;
		$nnih[$stap_nih][$ih][$ib-$off] = 1;
	    }
	}
	if ($ih >= 0 && $scaf_nih >= 0 && $ih != $scaf_nih) {
	    #print STDERR "scaf cross from H$ih to H$scaf_nih at $ib\n";
	    for my $off (0 .. 7) {
		$nnih[$ih][$scaf_nih][$ib+$off] = 1;
		$nnih[$scaf_nih][$ih][$ib+$off] = 1;
		$nnih[$ih][$scaf_nih][$ib-$off] = 1;
		$nnih[$scaf_nih][$ih][$ib-$off] = 1;
	    }
	}
    }
}

for my $ih (0 .. $#vstrands-1) {
    for my $jh ($ih+1 .. $#vstrands) {
	next unless defined $nnih[$ih][$jh];
	print STDERR "H$ih and H$jh are in contact\n";
	for my $ib (0 .. $#{$vstrands[$ih]->{stap}}) {
	    ## NOTE: $scaf_ri or $stap_ri can be undefined.
	    next unless defined $vstrands[$ih]{"scaf_ri"}[$ib][0];
	    next unless defined $vstrands[$ih]{"stap_ri"}[$ib][0];
	    next unless defined $vstrands[$jh]{"scaf_ri"}[$ib][0];
	    next unless defined $vstrands[$jh]{"stap_ri"}[$ib][0];

	    # skip if close to HJ.
	    next if defined $nnih[$ih][$jh][$ib];

	    my $scaf_ri = $vstrands[$ih]{"scaf_ri"}[$ib][0];
	    my $scaf_rj = $vstrands[$jh]{"scaf_ri"}[$ib][0];
	    my $stap_ri = $vstrands[$ih]{"stap_ri"}[$ib][0];
	    my $stap_rj = $vstrands[$jh]{"stap_ri"}[$ib][0];
	    printf FH "ri %d %d & a P\n", $scaf_ri, $stap_rj;
	    printf FH "ri %d %d & a P\n", $scaf_rj, $stap_ri;
	}
    }
}
print FH "q\n";
close FH;


####################################################
# chickenwire (CW)
open FH, ">chickenwire.for.make_ndx";
print FH "del 0 - 30\n";
print FH "!a C1' C2'  C3' C4' O4' C5'  P  O1P  O2P  O3'  O5' H*\n";
print FH "name 0 BASE\n";
# hbond pair
open FHB, ">hbonds.for.make_ndx";
print FHB "del 0 - 30\n";
print FHB "r DT* DC* & a N3\n";
print FHB "r DA* DG* & a N1\n";
print FHB "0 | 1\n";
print FHB "del 0\n";
print FHB "del 0\n";
print FHB "name 0 N1N3\n";

my $CWatom_num=0;
my @CWatoms;
my @CWres;
my @CWseg;
my %CWatoms_hash;
my @CWbonds;
my @CWjuncs;
my @CW_num_cross;   # crossover number at ib
for my $ih (0 .. $#vstrands) {
    my $dir = is_upstream($ih); # +1 for upstream, -1 for downstream.
    for my $ib (0 .. $#{$vstrands[$ih]->{stap}}) {
	next if ($vstrands[$ih]{stap}[$ib][1] == -1 
	      && $vstrands[$ih]{stap}[$ib][3] == -1);
	next if $vstrands[$ih]->{skip}[$ib] == -1;

	my ($stap_pih,$stap_pib,$stap_nih,$stap_nib) = @{$vstrands[$ih]{stap}[$ib]};
	my ($scaf_pih,$scaf_pib,$scaf_nih,$scaf_nib) = @{$vstrands[$ih]{scaf}[$ib]};

	# Atoms
	my @loops = (0 .. $vstrands[$ih]->{loop}[$ib]);
	@loops = reverse(@loops) if $dir == -1;
	my $inside_loop;
	for my $il (@loops) {
	    ++$CWatom_num;
	    ## NOTE: $scaf_ri or $stap_ri can be undefined.
	    my $scaf_ri = (defined $vstrands[$ih]{"scaf_ri"}[$ib][$il]) ? $vstrands[$ih]{"scaf_ri"}[$ib][$il] : "";
	    my $stap_ri = (defined $vstrands[$ih]{"stap_ri"}[$ib][$il]) ? $vstrands[$ih]{"stap_ri"}[$ib][$il] : "";

	    print FH "ri $scaf_ri $stap_ri & 0\n";
	    my $sname = "H$ih";
	    my $rname = "B$ib";
	    my $aname = "L$il";
	    push @CWatoms, $aname;
	    push @CWres, $rname;
	    push @CWseg, $sname;
	    my $key = "$sname$rname$aname";
	    $CWatoms_hash{$key} = $CWatom_num;
	    print FH "name $CWatom_num $key\n";

	    #################
	    # hbonds
	    if ((defined $vstrands[$ih]{"scaf_ri"}[$ib][$il]) && 
		(defined $vstrands[$ih]{"stap_ri"}[$ib][$il])) {
		my $x1 = $vstrands[$ih]{x1}[$ib][$il];
		my $y1 = $vstrands[$ih]{y1}[$ib][$il];
		my $z1 = $vstrands[$ih]{z1}[$ib][$il];
		my $x2 = $vstrands[$ih]{"x-1"}[$ib][$il];
		my $y2 = $vstrands[$ih]{"y-1"}[$ib][$il];
		my $z2 = $vstrands[$ih]{"z-1"}[$ib][$il];
		my $dx = $x1-$x2;
		my $dy = $y1-$y2;
		my $dz = $z1-$z2;
		my $dd = sqrt($dx*$dx+$dy*$dy+$dz*$dz);
		#printf STDERR "($x1,$y1,$z3)-($x2,$y2,$z3)=$dd\n";
		if ($dd > 4) {
		    printf STDERR "$scaf_ri-$stap_ri too far! ($x1,$y1)-($x2,$y2)\n";
		}
		print FHB "ri $scaf_ri & 0\n";
		print FHB "ri $stap_ri & 0\n";
		#print STDERR "hbonds: $scaf_ri $stap_ri\n";
	    }
	    #################

	    # bonds inside the loop
	    if (defined $inside_loop) {
		push @CWbonds, [$inside_loop,$key];
		$inside_loop = $key;
	    }
	}
	# ib-to-ib bonds.
	# only with next one to avoid double counting.
	if ($stap_nih != -1 && $stap_nib != -1) {
	    ($stap_nih,$stap_nib) = next_nonskip_ih_ib("stap",$stap_nih,$stap_nib);
	    push @CWbonds, ["H$ih"."B$ib"."L$loops[$#loops]","H$stap_nih"."B$stap_nib"."L0"];
	    #print STDERR "CW: stap cross H$ih B$ib L$loops[$#loops] to H$stap_nih B$stap_nib L0\n";
	    ++$CW_num_cross[$ib];
	}
	if ($scaf_nih != -1 && $scaf_nib != -1) {
	    ($scaf_nih,$scaf_nib) = next_nonskip_ih_ib("scaf",$scaf_nih,$scaf_nib);
	    push @CWbonds, ["H$ih"."B$ib"."L$loops[$#loops]","H$scaf_nih"."B$scaf_nib"."L0"];
	    #print STDERR "CW: stap cross H$ih B$ib L$loops[$#loops] to H$scaf_nih B$scaf_nib L0\n";
	    ++$CW_num_cross[$ib];
	}
    }
}
print FH "quit";
print FHB "quit";
close FH;
close FHB;
# write CW_num_cross
open FH, ">chickenwire.num.cross";
for my $i (0 .. $#CW_num_cross) {
    printf FH "%d %d\n", $i+1, (defined $CW_num_cross[$i])?$CW_num_cross[$i]:0;
}

## write psf
open FPSF, ">chickenwire.psf";
print FPSF "PSF EXT\n";
printf FPSF "%10d !NTITLE\n", 1;
printf FPSF "* Chickenwire representation made by cadnano2pdb.pl EXT\n\n";

printf FPSF "%10d !NATOM\n", $#CWatoms+1;
#fmt02='(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6)'

my $index=0;
my $resid=1;
for my $a (0 .. $#CWatoms) {
    printf FPSF "%10d ", ++$index;
    printf FPSF "%-8s ", $CWseg[$a];
    printf FPSF "%-8d ", $index;
    printf FPSF "%-8s ", $CWres[$a];
    printf FPSF "%-6s ", $CWatoms[$a];
    printf FPSF "%6d " , 6;
    print FPSF  "  0.000000       1.00800";
    print FPSF  "           0   0.00000     -0.301140E-02\n";
}

## replace CWbonds using hash
my @CWbonds_num;
for my $i (0 .. $#CWbonds) {
    #print STDERR "key=$CWbonds[$i][0] $CWatoms_hash{$CWbonds[$i][0]}\n";
    #print STDERR "key=$CWbonds[$i][1] $CWatoms_hash{$CWbonds[$i][1]}\n";
    if (!defined $CWatoms_hash{$CWbonds[$i][0]}) {
	print STDERR "Undefined $CWbonds[$i][0]\n";
	next;
    }
    if (!defined $CWatoms_hash{$CWbonds[$i][1]}) {
	print STDERR "Undefined $CWbonds[$i][1]\n";
	next;
    }
    my $bi = $CWatoms_hash{$CWbonds[$i][0]};
    my $bj = $CWatoms_hash{$CWbonds[$i][1]};
    printf STDERR "Bond($bi,$bj) out of range\n" if ($bi < 1 || $bi > ($#CWatoms+1));
    printf STDERR "Bond($bi,$bj) out of range\n" if ($bj < 1 || $bj > ($#CWatoms+1));
    push @CWbonds_num, ($bi,$bj);
}
printf FPSF "\n%10d !NBOND: bonds\n", ($#CWbonds_num+1)/2;
while ($#CWbonds_num >= 7) {
    #printf STDERR "%10d%10d%10d%10d%10d%10d%10d%10d\n", splice(@CWbonds_num, 0, 8);
    printf FPSF "%10d%10d%10d%10d%10d%10d%10d%10d\n", splice(@CWbonds_num, 0, 8);
}
while ($#CWbonds_num >= 0) {
    printf FPSF "%10d", shift(@CWbonds_num);
}
print FPSF "\n";

printf FPSF "\n%10d !NTHETA: angles\n", 0;
printf FPSF "\n%10d !NPHI: dihedrals\n", 0;
printf FPSF "\n%10d !NIMPHI: impropers\n", 0;
printf FPSF "\n%10d !NDON: donors\n", 0;
printf FPSF "\n%10d !NACC: acceptors\n", 0;
close FPSF;

exit;

####################################################
# H-bonds
#open FHGMX, ">hbonds.itp";
#open FHNAMD, ">hbonds.extrabonds";
#print FHGMX "[ bonds ]\n";
#for my $ih (0 .. $#vstrands) {
#    my $c1p_stap_old;
#    my $c1p_scaf_old;
#    for my $ib (0 .. $#{$vstrands[$ih]->{stap}}) {
#	#printf STDERR "ih=%d ib=%d [%d %d %d %d]\n", $ih,$ib,@{$vstrands[$ih]{stap}[$ib]};
#	next if ($vstrands[$ih]{stap}[$ib][1] == -1 && $vstrands[$ih]{stap}[$ib][3] == -1);
#	next if $vstrands[$ih]->{skip}[$ib] == -1;
#
#	for my $il (0 .. $vstrands[$ih]->{loop}[$ib]) {
#	    #### H-bonds
#	    printf FHGMX "%8d %8d %8d %5.2f %5.2f %5.2f %5.2f ; ih=%d ib=%d major\n",
#		    @{$vstrands[$ih]{hb_major}[$ib][$il]},10,0.0,0.30,9.0,1000,$ih,$ib;
#	    printf FHGMX "%8d %8d %8d %5.2f %5.2f %5.2f %5.2f ; ih=%d ib=%d center\n",
#		    @{$vstrands[$ih]{hb_center}[$ib][$il]},10,0.0,0.30,9.0,1000,$ih,$ib;
#	    printf FHGMX "%8d %8d %8d %5.2f %5.2f %5.2f %5.2f ; ih=%d ib=%d minor\n",
#		    @{$vstrands[$ih]{hb_minor}[$ib][$il]},10,0.0,0.30,9.0,1000,$ih,$ib;
#	    printf FHNAMD "bond%10d%10d%10.3g%10.3g\n",
#		    $vstrands[$ih]{hb_major}[$ib][$il][0]-1,
#		    $vstrands[$ih]{hb_major}[$ib][$il][1]-1, 1, 2.8;
#	    printf FHNAMD "bond%10d%10d%10.3g%10.3g\n",
#		    $vstrands[$ih]{hb_center}[$ib][$il][0]-1,
#		    $vstrands[$ih]{hb_center}[$ib][$il][1]-1, 1, 2.8;
#	    printf FHNAMD "bond%10d%10d%10.3g%10.3g\n",
#		    $vstrands[$ih]{hb_minor}[$ib][$il][0]-1,
#		    $vstrands[$ih]{hb_minor}[$ib][$il][1]-1, 1, 2.8;
#
#	    #### C1'
#	    if (defined $c1p_stap_old) {
#		printf FHGMX "%8d %8d %8d %5.2f %5.2f %5.2f %5.2f ; C1' staple\n",
#		    $c1p_stap_old,$vstrands[$ih]{c1p_stap}[$ib][$il],
#		    10,0.0,0.60,9.0,1000;
#		printf FHGMX "%8d %8d %8d %5.2f %5.2f %5.2f %5.2f ; C1' scafle\n",
#		    $c1p_scaf_old,$vstrands[$ih]{c1p_scaf}[$ib][$il],
#		    10,0.0,0.60,9.0,1000;
#	    }
#	    $c1p_stap_old = $vstrands[$ih]{c1p_stap}[$ib][$il];
#	    $c1p_scaf_old = $vstrands[$ih]{c1p_scaf}[$ib][$il];
#	}
#    }
#}
#close FHGMX;
#close FHNAMD;


####################################################
# H-bonds
# Delete hbonds in print_pdb.
open FH, ">hbonds.itp";
open FHNAMD, ">hbonds.extrabonds";
print FH "[ bonds ]\n";
for my $ih (0 .. $#vstrands) {
    for my $ib (0 .. $#{$vstrands[$ih]->{stap}}) {
	next if ($vstrands[$ih]{stap}[$ib][1] == -1 
	      && $vstrands[$ih]{stap}[$ib][3] == -1);
	next if $vstrands[$ih]->{skip}[$ib] == -1;

	for my $il (0 .. $vstrands[$ih]->{loop}[$ib]) {
	    my $stap_nt = $vstrands[$ih]{"stap_nt"}[$ib][$il];
	    my $a1;
	    my $a2;
	    if ($stap_nt =~ "(A|G)") {
		$a1 = ' N1 ';
		$a2 = ' N3 ';
	    }
	    else {
		$a1 = ' N3 ';
		$a2 = ' N1 ';
	    }
	    my $seg = int($ib/$nblock);
	    printf FHNAMD "%8d %8d %s %s\n", 
				$vstrands[$ih]{"stap_atoms"}[$ib][$il]{$a1}-1,
				$vstrands[$ih]{"scaf_atoms"}[$ib][$il]{$a2}-1,
				"H$ih"."S$seg", "$a1-$a2";
	    ## atom 
	    #push @buf, $vstrands[$ih]{"scaf$aa"}[$ib][$il]
	    #    if defined $vstrands[$ih]{"scaf$aa"}[$ib][$il];
	    #
	    #push @buf, $vstrands[$ih]{"stap$aa"}[$ib][$il]
	    #    if defined $vstrands[$ih]{"stap$aa"}[$ib][$il];
	}
    }
}
close FH;
close FHNAMD;

####################################################
# Holliday junction.
# Print residue indices that define four vectors (from junction to arm) for make_ndx.
# cat junction.for.make_ndx | make_ndx -f tpr -o junction.ndx
#            8    5
#            ^    |
#            |    V
#  (ih,ib+1) 7 <- 6 (nih,ib+1)
#  (ih,ib)   2 -> 3 (nih,ib)
#            ^    |
#            |    V
#            1    4
open FH, ">junction.for.make_ndx";
print FH "del 0-30\n";
for my $ih (0 .. $#vstrands) {
    next if is_upstream($ih) == -1;
    my @HJ_queue=();
    my $dir;	# +1 for upstream, -1 for downstream.
    for my $ib (0 .. $#{$vstrands[$ih]->{stap}}) {
	next if ($vstrands[$ih]{stap}[$ib][1] == -1 
	      && $vstrands[$ih]{stap}[$ib][3] == -1);
	next if $vstrands[$ih]->{skip}[$ib] == -1;

	## We only count upstream helix to avoid double counting HJ.
	my $pih = $vstrands[$ih]{stap}[$ib][0];
	my $pib = $vstrands[$ih]{stap}[$ib][1];
	my $nih = $vstrands[$ih]{stap}[$ib][2];
	my $nib = $vstrands[$ih]{stap}[$ib][3];

	#################################
	# CAVEAT: loop NTs are ordered OK only for upstream:il=0,1,2,... 
	# Order of loop NTs must be reverse for downstream:il=...,1,0.
	#for my $il (0 .. $vstrands[$ih]->{loop}[$ib]) {
	#    # push all loop NTs because crossing is impossible in a loop.
	#    push @HJ_queue, [$ih,$ib,$il];
	#}

	## check if you're in a junction.
	next if $nih == -1;
	next if $pih == -1;
	next if ($nih == $ih);

	#################################
	# Found a junction.
	# 2 = $ih,$ib,last_il
	# 3 = $nih,$ib,0
	# 7 = $ih,$ib+1,0
	# 6 = $nih,$ib+1,last_il
	print STDERR "HJ ($ih,$ib)-($nih,$nib)\n";
	my @a2to1 = trace_helix_downstream($ih,$ib);
	my @a3to4 = trace_helix_downstream($nih,$ib);
	my @a6to5 = trace_helix_upstream($nih,$ib+1);
	my @a7to8 = trace_helix_upstream($ih,$ib+1);
	my $r1a = $vstrands[$a2to1[3][0]]->{stap_ri}[$a2to1[3][1]][$a2to1[3][2]];
	my $r1b = $vstrands[$a2to1[3][0]]->{scaf_ri}[$a2to1[3][1]][$a2to1[3][2]];
	my $r2a = $vstrands[$a2to1[0][0]]->{stap_ri}[$a2to1[0][1]][$a2to1[0][2]];
	my $r2b = $vstrands[$a2to1[0][0]]->{scaf_ri}[$a2to1[0][1]][$a2to1[0][2]];
	my $r3a = $vstrands[$a3to4[0][0]]->{stap_ri}[$a3to4[0][1]][$a3to4[0][2]];
	my $r3b = $vstrands[$a3to4[0][0]]->{scaf_ri}[$a3to4[0][1]][$a3to4[0][2]];
	my $r4a = $vstrands[$a3to4[3][0]]->{stap_ri}[$a3to4[3][1]][$a3to4[3][2]];
	my $r4b = $vstrands[$a3to4[3][0]]->{scaf_ri}[$a3to4[3][1]][$a3to4[3][2]];
	my $r5a = $vstrands[$a6to5[3][0]]->{stap_ri}[$a6to5[3][1]][$a6to5[3][2]];
	my $r5b = $vstrands[$a6to5[3][0]]->{scaf_ri}[$a6to5[3][1]][$a6to5[3][2]];
	my $r6a = $vstrands[$a6to5[0][0]]->{stap_ri}[$a6to5[0][1]][$a6to5[0][2]];
	my $r6b = $vstrands[$a6to5[0][0]]->{scaf_ri}[$a6to5[0][1]][$a6to5[0][2]];
	my $r7a = $vstrands[$a7to8[0][0]]->{stap_ri}[$a7to8[0][1]][$a7to8[0][2]];
	my $r7b = $vstrands[$a7to8[0][0]]->{scaf_ri}[$a7to8[0][1]][$a7to8[0][2]];
	my $r8a = $vstrands[$a7to8[3][0]]->{stap_ri}[$a7to8[3][1]][$a7to8[3][2]];
	my $r8b = $vstrands[$a7to8[3][0]]->{scaf_ri}[$a7to8[3][1]][$a7to8[3][2]];
	### 
	# Reorder ri to be consistent with indices in g_traj2holliday.pl
	# 1 2 4 3 5 6 8 7
	print FH "ri $r1a $r1b\n";
	print FH "ri $r2a $r2b\n";
	print FH "ri $r4a $r4b\n";
	print FH "ri $r3a $r3b\n";
	print FH "ri $r5a $r5b\n";
	print FH "ri $r6a $r6b\n";
	print FH "ri $r8a $r8b\n";
	print FH "ri $r7a $r7b\n";
    }
}
print FH "q\n";
close FH;

exit;

####################################################
# ENM bonds
open FH, ">enm_bdna.itp";
open FHNAMD, ">enm_bdna.namd.extrabonds";
print FH "[ bonds ]\n";
for my $ih (0 .. $#vstrands) {
    #### determine directionn of scaffold.
    # For example, [0,5,0,7] means scaffold moves from bottom to top (+z direction).
    # staple flows in opposite direction.
    my $dir;
    for my $ib (0 .. $#{$vstrands[$ih]->{scaf}}) {
	my $five = $vstrands[$ih]{scaf}[$ib][1];
	my $three = $vstrands[$ih]{scaf}[$ib][3];
	if ($five != -1 && $three != -1) {
	    $dir = ($five > $three) ? 1 : -1;
	    last;
	}
    }

    my $old_ib;
    my $old_il;
    for my $ib (0 .. $#{$vstrands[$ih]->{stap}}) {
	next if ($vstrands[$ih]{stap}[$ib][1] == -1 
	      && $vstrands[$ih]{stap}[$ib][3] == -1);
	next if $vstrands[$ih]->{skip}[$ib] == -1;

	for my $il (0 .. $vstrands[$ih]->{loop}[$ib]) {
	    if (defined $old_ib) { #### print ENM bonds
		if ($dir > 0) {
		    print_enm_stack($ih,$old_ib,$old_il,$ib,$il);
		}
		else {
		    print_enm_stack($ih,$ib,$il,$old_ib,$old_il);
		}
	    }
	    print_enm_basepair($ih,$ib,$il);

	    $old_ib = $ib;
	    $old_il = $il;
	}
    }
}
close FH;
close FHNAMD;

####################################################
# Index for base atoms
# $nblock = 7 (honeycomb) or 8 (square)
open FH, ">segments.ndx";
open FHRES, ">segments_residue.ndx";
my $fhres_nn=1;
print FHRES "del 0 - 30\n";
print FHRES "!a C1' C2'  C3' C4' O4' C5'  P  O1P  O2P  O3'  O5' H*\n";
print FHRES "name 0 BASE\n";
my $CGatom_num = -1;
my @CGatoms;
my @CGbonds;
for my $ih (0 .. $#vstrands) {
    my @buf;
    my @buf_ri;
    my $seg=0;
    for my $ib (0 .. $#{$vstrands[$ih]->{stap}}) {
	## flush if any
	if (($ib % $nblock) == 0 && @buf ) {
	    # atom
	    my $nn=0;
	    print FH "[ H$ih"."S$seg ]\n";
	    for my $i (sort { $a <=> $b } @buf) {
		printf FH "%8d ", $i;
		print FH "\n" if (++$nn % 10) == 0;
	    }
	    print FH "\n" unless ($nn % 10) == 0;
	    @buf = ();
	    print FH "\n";

	    # residue
	    print FHRES "ri ";
	    for my $i (sort { $a <=> $b } @buf_ri) {
		printf FHRES "%d ", $i;
	    }
	    print FHRES "& \"BASE\"\n";
	    print FHRES "name $fhres_nn H$ih"."S$seg\n";
	    ++$fhres_nn;
	    @buf_ri = ();

	    ## PSF
	    push @CGatoms, "H$ih"."S$seg";
	    ++$seg;
	}

	next if ($vstrands[$ih]{stap}[$ib][1] == -1 
	      && $vstrands[$ih]{stap}[$ib][3] == -1);
	next if $vstrands[$ih]->{skip}[$ib] == -1;

	for my $il (0 .. $vstrands[$ih]->{loop}[$ib]) {
	    for my $aa (@baseatoms) {
		## atom 
		push @buf, $vstrands[$ih]{"scaf$aa"}[$ib][$il]
		    if defined $vstrands[$ih]{"scaf$aa"}[$ib][$il];

		push @buf, $vstrands[$ih]{"stap$aa"}[$ib][$il]
		    if defined $vstrands[$ih]{"stap$aa"}[$ib][$il];
	    }
	    ## residue
	    push @buf_ri, $vstrands[$ih]{"scaf_ri"}[$ib][$il]
		if defined $vstrands[$ih]{"scaf_ri"}[$ib][$il];

	    push @buf_ri, $vstrands[$ih]{"stap_ri"}[$ib][$il]
		if defined $vstrands[$ih]{"stap_ri"}[$ib][$il];
	}
    }
    ## flush if any
    if ( @buf ) {
	# atom
	my $nn=0;
	print FH "[ H$ih"."S$seg ]\n";
	for my $i (sort { $a <=> $b } @buf) {
	    printf FH "%8d ", $i;
	    print FH "\n" if (++$nn % 10) == 0;
	}
	print FH "\n" unless ($nn % 10) == 0;
	@buf = ();
	print FH "\n";

	# residue
	print FHRES "ri ";
	for my $i (sort { $a <=> $b } @buf_ri) {
	    printf FHRES "%d ", $i;
	}
	print FHRES "& \"BASE\"\n";
	print FHRES "name $fhres_nn H$ih"."S$seg\n";
	++$fhres_nn;
	@buf_ri = ();

	## PSF
	push @CGatoms, "H$ih"."S$seg";
    }

    ## PSF
    for my $i (($CGatom_num+1) .. ($#CGatoms-1)) {
	push @CGbonds, $i+1, $i+2;
    }
    $CGatom_num = $#CGatoms;
}
close FH;
print FHRES "q\n";
close FHRES;
## write psf
open FPSF, ">segments.psf";
printf FPSF "%10d !NTITLE\n", 1;
printf FPSF "* topol.tpr converted to PSF EXT\n\n";

printf FPSF "%10d !NATOM\n", $#CGatoms+1;
#fmt02='(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6)'

$index=0;
$resid=1;
for my $a (0 .. $#CGatoms) {
    printf FPSF "%10d ", ++$index;
    printf FPSF "%-8s ", $CGatoms[$a];
    printf FPSF "%-8d ", $index;
    printf FPSF "%-8s ", $CGatoms[$a];
    printf FPSF "%-6s ", $CGatoms[$a];
    printf FPSF "%6d " , 6;
    print FPSF  "  0.000000       1.00800";
    print FPSF  "           0   0.00000     -0.301140E-02\n";
}

printf FPSF "\n%10d !NBOND: bonds\n", ($#CGbonds+1)/2;
while ($#CGbonds >= 7) {
    printf FPSF "%10d%10d%10d%10d%10d%10d%10d%10d\n", splice(@CGbonds, 0, 8);
}
while ($#CGbonds >= 0) {
    printf FPSF "%10d", shift(@CGbonds);
}
print FPSF "\n";

printf FPSF "\n%10d !NTHETA: angles\n", 0;
printf FPSF "\n%10d !NPHI: dihedrals\n", 0;
printf FPSF "\n%10d !NIMPHI: impropers\n", 0;
printf FPSF "\n%10d !NDON: donors\n", 0;
printf FPSF "\n%10d !NACC: acceptors\n", 0;
close FPSF;

exit;

###############################################
# subroutines 
###############################################

sub print_enm_stack {
    my $ih  = shift;
    my $ib5 = shift;
    my $il5 = shift;
    my $ib3 = shift;
    my $il3 = shift;

    ############################
    ### scaffold
    ############################
    my $nt5 = wc($vstrands[$ih]{stap_nt}[$ib5][$il5]);
    my $nt3 = wc($vstrands[$ih]{stap_nt}[$ib3][$il3]);
    my $key = "D".$nt5."D".$nt3;
    $key =~ tr/GC/AT/;
    #print STDERR "$key ($ih),($ib5,$il5,$nt5),($ib3,$il3,$nt3)\n";

    for my $a5 (@baseatoms) {
	next unless defined $vstrands[$ih]{"scaf$a5"}[$ib5][$il5];
	for my $a3 (@baseatoms) {
	    next unless defined $vstrands[$ih]{"scaf$a3"}[$ib3][$il3];
	    #print STDERR "$nt5($a5) $nt3($a3)\n";
	    printf FH "%10d%10d%10d%10.3g%10d%s\n",
		    $vstrands[$ih]{"scaf$a5"}[$ib5][$il5],
		    $vstrands[$ih]{"scaf$a3"}[$ib3][$il3],
		    6, $enm_stack{$key}{$a5}{$a3}*0.1, 100,
		    " ; $nt5($a5) - $nt3($a3) stack";
	    ## NAMD extrabonds
	    printf FHNAMD "bond%10d%10d%10.3g%10.3g\n",
		    $vstrands[$ih]{"scaf$a5"}[$ib5][$il5]-1,
		    $vstrands[$ih]{"scaf$a3"}[$ib3][$il3]-1,
		    0.2, $enm_stack{$key}{$a5}{$a3};
	}
    }
    ############################
    ### staple
    ############################
    ($ib5,$ib3) = ($ib3,$ib5);
    ($il5,$il3) = ($il3,$il5);
    $nt5 = $vstrands[$ih]{stap_nt}[$ib5][$il5];
    $nt3 = $vstrands[$ih]{stap_nt}[$ib3][$il3];
    $key = "D".$nt5."D".$nt3;
    $key =~ tr/GC/AT/;
    #print STDERR "$key ($ih),($ib5,$il5,$nt5),($ib3,$il3,$nt3)\n";

    for my $a5 (@baseatoms) {
	next unless defined $vstrands[$ih]{"stap$a5"}[$ib5][$il5];
	for my $a3 (@baseatoms) {
	    next unless defined $vstrands[$ih]{"stap$a3"}[$ib3][$il3];
	    #print STDERR "$nt5($a5) $nt3($a3)\n";
	    printf FH "%10d%10d%10d%10.3g%10d%s\n",
		    $vstrands[$ih]{"stap$a5"}[$ib5][$il5],
		    $vstrands[$ih]{"stap$a3"}[$ib3][$il3],
		    6, $enm_stack{$key}{$a5}{$a3}*0.1, 100,
		    " ; $nt5($a5) - $nt3($a3) stack";
	    ## NAMD extrabonds
	    printf FHNAMD "bond%10d%10d%10.3g%10.3g\n",
		    $vstrands[$ih]{"stap$a5"}[$ib5][$il5]-1,
		    $vstrands[$ih]{"stap$a3"}[$ib3][$il3]-1,
		    0.2, $enm_stack{$key}{$a5}{$a3};
	}
    }
}

sub print_enm_basepair {
    my $ih = shift;
    my $ib = shift;
    my $il = shift;

    my $nt2 = $vstrands[$ih]{stap_nt}[$ib][$il];
    my $nt1 = wc($nt2);
    my $key = "D".$nt1."D".$nt2;
    $key =~ tr/GC/AT/;

    for my $a1 (@baseatoms) {
	next unless defined $vstrands[$ih]{"scaf$a1"}[$ib][$il];
	for my $a2 (@baseatoms) {
	    next unless defined $vstrands[$ih]{"stap$a2"}[$ib][$il];
	    #print STDERR "$nt5($a1) $nt3($a2)\n";
	    printf FH "%10d%10d%10d%10.3g%10d%s\n",
		    $vstrands[$ih]{"scaf$a1"}[$ib][$il],
		    $vstrands[$ih]{"stap$a2"}[$ib][$il],
		    6, $enm_bp{$key}{$a1}{$a2}*0.1, 500,
		    " ; $nt1($a1) - $nt2($a2) bp";
	}
    }
}

sub assign_seq_staple {
    my ($ih,$ib) = @_;

    # find the corresponding sequence in csv.
    my $pseq;	# pointer to sequence array.
    for my $i (0 .. $#csv) {
	next unless ($csv[$i]{'5helix'} == $ih);
	next unless ($csv[$i]{'5base'} == $ib);

	# now we found it.
	$pseq = $csv[$i]{seq};
	$nseg = $i;
	last;
    }
    die "Can't find a sequence for helix=$ih,base=$ib" unless defined $pseq;

    ## trace the staple.
    do {
	my $dir = get_direction_stap($ih);	# upstream=1, downstream=-1
	# next base
	#print STDERR "ih=$ih,ib=$ib\n";
	my ($pih,$pib,$nih,$nib) = @{$vstrands[$ih]{stap}[$ib]};

	if ($vstrands[$ih]{skip}[$ib] == -1) {	# skip
	    # Add sequence to vstrands.
	    #$vstrands[$ih]{stap_nt}[$ib] = 'X';
	}
	else {
	    my @list_il = (0 .. $vstrands[$ih]{loop}[$ib]); # upstream
	    @list_il = reverse(@list_il) if ($dir < 0);	# downstream
	    for my $il (@list_il) {
		my $nt = shift @{$pseq};
		# Add sequence to vstrands.
		$vstrands[$ih]{stap_nt}[$ib][$il] = $nt;
		$vstrands[$ih]{scaf_nt}[$ib][$il] = wc($nt);
	    }
	}

	# next base
	($ih,$ib) = ($nih,$nib);
    } until ($ib == -1);
}

sub check_seq_scaf {
    my ($ih,$ib) = @_;

    ## trace the scaffold.
    do {
	my $dir = get_direction_scaf($ih);	# upstream=1, downstream=-1
	# next base
	my ($pih,$pib,$nih,$nib) = @{$vstrands[$ih]{scaf}[$ib]};

	if ($vstrands[$ih]{skip}[$ib] == -1) {	# skip
	}
	else {
	    my @list_il = (0 .. $vstrands[$ih]{loop}[$ib]); # upstream
	    @list_il = reverse(@list_il) if ($dir < 0);	# downstream
	    for my $il (@list_il) {
		unless (defined $vstrands[$ih]{scaf_nt}[$ib][$il]) {
		    die "Undefined scaffold sequence at ih=$ih,ib=$ib,il=$il\n";
		}
	    }
	}
	# next base
	($ih,$ib) = ($nih,$nib);
    } until ($ib == -1);
}

sub check_seq_stap {
    my ($ih,$ib) = @_;

    ## trace the staple.
    do {
	my $dir = get_direction_stap($ih);	# upstream=1, downstream=-1
	# next base
	my ($pih,$pib,$nih,$nib) = @{$vstrands[$ih]{stap}[$ib]};

	if ($vstrands[$ih]{skip}[$ib] == -1) {	# skip
	}
	else {
	    my @list_il = (0 .. $vstrands[$ih]{loop}[$ib]); # upstream
	    @list_il = reverse(@list_il) if ($dir < 0);	# downstream
	    for my $il (@list_il) {
		unless (defined $vstrands[$ih]{stap_nt}[$ib][$il]) {
		    #die "Undefined staple sequence at ih=$ih,ib=$ib,il=$il\n";
		    print STDERR "Using T for undefined staple sequence at ih=$ih,ib=$ib,il=$il\n";
		    $vstrands[$ih]{stap_nt}[$ib][$il] = 'T';
		}
	    }
	}
	# next base
	($ih,$ib) = ($nih,$nib);
    } until ($ib == -1);
}

sub trace_staple {
    my ($ih,$ib) = @_;

    ## trace the strand.
    my $resi=1;
    do {
	my $dir = get_direction_stap($ih);	# upstream=1, downstream=-1
	# next base
	#print STDERR "ih=$ih,ib=$ib\n";
	my ($pih,$pib,$nih,$nib) = @{$vstrands[$ih]{stap}[$ib]};

	if ($vstrands[$ih]{skip}[$ib] == -1) {	# skip
	    # Add sequence to vstrands.
	    $vstrands[$ih]{stap_nt}[$ib] = 'X';
	}
	else {
	    # coordinates
	    my $hrow = $vstrands[$ih]{row};	# helix row
	    my $hcol = $vstrands[$ih]{col};	# helix column

	    my @list_il = (0 .. $vstrands[$ih]{loop}[$ib]); # upstream
	    @list_il = reverse(@list_il) if ($dir < 0);	# downstream
	    for my $il (@list_il) {
		# Add sequence to vstrands.
		#push @{$vstrands[$ih]{stap_nt}[$ib]}, $nt;
		my $nt = $vstrands[$ih]{stap_nt}[$ib][$il];
		# global residue count
		$vstrands[$ih]{stap_ri}[$ib][$il] = ++$global_ri;
		#print STDERR "stap: $ih $ib $il $global_ri\n";

		my $theta = $vstrands[$ih]{theta}[$ib][$il];	# 
		my $z = $vstrands[$ih]{z}[$ib][$il];	# 

		my $resn = na1to3($nt);
		if ($resi == 1) {
		    $resn = $resn."5";	# 5-ter
		}
		elsif ($nib == -1) {
		    $resn = $resn."3";	# 3-ter
		}
		my $seg = sprintf "P%03d", $nseg;

		my %atoms = print_pdb($hrow,$hcol,$ih,$ib,$il,$nt,$dir,$seg,$resi,$resn,$theta,$z);
		$vstrands[$ih]{stap_atoms}[$ib][$il] = \%atoms;
		#### H-bond atom indices
		#push @{$vstrands[$ih]{hb_major}[$ib][$il]}, $enm_atoms[0];	# 
		#push @{$vstrands[$ih]{hb_center}[$ib][$il]}, $enm_atoms[1];	# 
		#push @{$vstrands[$ih]{hb_minor}[$ib][$il]}, $enm_atoms[2];	# 
		#$vstrands[$ih]{c1p_stap}[$ib][$il] = $enm_atoms[3];	# 
		### base atoms
		#for my $ba (0 .. $#baseatoms) {
		#    if (defined $enm_atoms[4+$ba]) {
		#	$vstrands[$ih]{"stap$baseatoms[$ba]"}[$ib][$il] = $enm_atoms[4+$ba];
		#    }
		#}

		++$resi;
	    }
	}

	# next base
	($ih,$ib) = ($nih,$nib);
    } until ($ib == -1);
    #print "TER\n";
    ++$ch;
    #$ch -= 26 if ($ch == 91);
}

sub trace_scaffold {
    my ($ih,$ib) = @_;

    my $resi=0;
    ## trace the strand.
    do {
	# next base
	my $dir = get_direction_scaf($ih);	# upstream=1, downstream=-1
	my ($pih,$pib,$nih,$nib) = @{$vstrands[$ih]{scaf}[$ib]};
	#printf STDERR "[%d,%d,%d,%d]\n", $pih,$pib,$nih,$nib;

	if ($vstrands[$ih]{skip}[$ib] == -1) {	# skip
	}
	else {
	    # coordinates
	    my $hrow = $vstrands[$ih]{row};	# helix row
	    my $hcol = $vstrands[$ih]{col};	# helix column

	    my @list_il = (0 .. $vstrands[$ih]{loop}[$ib]);
	    @list_il = reverse(@list_il) if ($dir < 0);
	    for my $il (@list_il) {
		++$resi;

		my $nt = $vstrands[$ih]{scaf_nt}[$ib][$il];
		# global residue count
		$vstrands[$ih]{scaf_ri}[$ib][$il] = ++$global_ri;
		#print STDERR "scaf: $ih $ib $il $global_ri\n";

		my $theta = $vstrands[$ih]{theta}[$ib][$il];	# 
		my $z = $vstrands[$ih]{z}[$ib][$il];	# 

		my $resn = na1to3($nt);
		if ($resi == 1) {
		    $resn = $resn."5";	# 5-ter
		}
		elsif ($nib == -1) {
		    $resn = $resn."3";	# 3-ter
		}

		my %atoms = print_pdb($hrow,$hcol,$ih,$ib,$il,$nt,$dir,"SCAF",$resi,$resn,$theta,$z);
		$vstrands[$ih]{scaf_atoms}[$ib][$il] = \%atoms;
		#### H-bond atom indices
		#push @{$vstrands[$ih]{hb_major}[$ib][$i]}, $enm_atoms[0];	# 
		#push @{$vstrands[$ih]{hb_center}[$ib][$i]}, $enm_atoms[1];	# 
		#push @{$vstrands[$ih]{hb_minor}[$ib][$i]}, $enm_atoms[2];	# 
		#$vstrands[$ih]{c1p_scaf}[$ib][$il] = $enm_atoms[3];	# 
		#for my $ba (0 .. $#baseatoms) {
		#    if (defined $enm_atoms[4+$ba]) {
		#	$vstrands[$ih]{"scaf$baseatoms[$ba]"}[$ib][$il] = $enm_atoms[4+$ba];
		#    }
		#}
	    }
	}

	# next base
	($ih,$ib) = ($nih,$nib);
    } until ($ib == -1);
    print "TER\n";
    ++$ch;
    #$ch -= 26 if ($ch == 91);

    # update $scaf_resi for the next scaffold strand.
    $scaf_resi += $resi;
}

##########################################################
# In case of multiple scaffolds, we assign sequence this way:
# 1. concatenate scaffold strands in order of 5' end apperance starting from ih=0,ir=0.
# 2. assign sequence to the concatenated scaffolds in 1.
# scaf_resi is the global variable for this.
sub assign_seq_scaffold {
    my ($ih,$ib) = @_;

    my $resi=0;
    ## trace the strand.
    do {
	# next base
	my $dir = get_direction_scaf($ih);	# upstream=1, downstream=-1
	my ($pih,$pib,$nih,$nib) = @{$vstrands[$ih]{scaf}[$ib]};
	#printf STDERR "[%d,%d,%d,%d]\n", $pih,$pib,$nih,$nib;

	if ($vstrands[$ih]{skip}[$ib] == -1) {	# skip
	}
	else {
	    my @list_il = (0 .. $vstrands[$ih]{loop}[$ib]);
	    @list_il = reverse(@list_il) if ($dir < 0);
	    for my $il (@list_il) {
		++$resi;
		if (($scaf_resi+$resi-1) > $#{$scafseq{$scafseqname}}) {
		    printf STDERR "Error: %d is more than $scafseqname seq\n", ($scaf_resi+$resi-1);
		    die;
		}
		my $nt = $scafseq{$scafseqname}[$scaf_resi+$resi-1];
		$vstrands[$ih]{scaf_nt}[$ib][$il] = $nt;
	    }
	}

	# next base
	($ih,$ib) = ($nih,$nib);
    } until ($ib == -1);

    # update $scaf_resi for the next scaffold strand.
    $scaf_resi += $resi;
}

sub print_pdb {
    my $ix = shift; # row
    my $iy = shift; # col
    my $ih = shift;
    my $ib = shift;
    my $il = shift;
    my $seq = shift; # ATGC
    my $dir = shift; # upstream (1) or downstream (-1)
    my $seg = shift;	# segment id
    my $resid = shift;
    my $resn = shift;
    my $theta = shift;
    my $z = shift;

    #### CHARMM 3TER rearrange atoms... dirty.
    #### Always doublecheck atom indices made by this script & the PSF file!!
    if ($forcefield eq "charmm" && $resn =~ "3") {
	$seq .= "3";
    }
    my %atoms;
    for my $i (@{$pdb{$seq}}) {
	my $aname = $i->[0];
	my @xyz = @{$i->[1]};

	# 5' terminal
	if ($resid == 1) {
	    if ($aname =~ / P/) {
		if ($forcefield eq "amber") {
		    $aname = " H5T";
		}
		else {
		    ## in case of CHARMM, put H5T before O5'. Dirty...
		    next;
		}
	    }
	    elsif ($aname =~ /O[12]P/) {
		next;
	    }
	}

	# flip PDB if downstream.
	@xyz = flip_x(@xyz) if ($dir == -1);

	if ($lattice eq "honeycomb") {
	    # rotate wrt z
	    #@xyz = rotate_z(@xyz,$ib*360*2/21);
	    @xyz = rotate_z(@xyz,$theta);
	    ########################################################
	    # 180 +- ??. Fine tuning needed.
	    @xyz = rotate_z(@xyz,170);

	    # translate helix
	    if ((($iy+$ix)%2) == 0) {
		$xyz[0] -= $ix*(1+1/sqrt(3))*$dd;
	    }
	    else {
		$xyz[0] -= ($ix+($ix+1)/sqrt(3))*$dd;
	    }
	    $xyz[1] += $iy*sqrt(3)/2.0*$dd;
	    #$xyz[2] += $ib*3.4;
	    $xyz[2] += $z;
	}
	elsif ($lattice eq "square") {
	    # rotate wrt z
	    #@xyz = rotate_z(@xyz,$ib*360*3/32);
	    @xyz = rotate_z(@xyz,$theta);
	    ########################################################
	    # 180 +- ??. Fine tuning needed.
	    @xyz = rotate_z(@xyz,15);

	    # translate helix
	    $xyz[0] -= $ix*$dd;
	    $xyz[1] += $iy*$dd;
	    #$xyz[2] += $ib*3.4;
	    $xyz[2] += $z;
	}

	# CHARMM H5T before O5'
	if ($forcefield eq "charmm" && $resn =~ /D.5/ && $aname =~ /O5'/) {
	    ++$global_aid;
	    $atoms{$aname} = $global_aid;
	    printf "%-6s%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n",
    		"ATOM",0, " H5T", $resn, num2chr($ch), $resid, 
		$xyz[0]+rand, $xyz[1]+rand, $xyz[2]+rand,1,0, $seg;
	}

	++$global_aid;
	$atoms{$aname} = $global_aid;
        printf "%-6s%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n",
    		"ATOM",0, $aname, $resn, num2chr($ch), $resid, @xyz, 1,0, $seg;

	# N1|N3 atom position
	if ( ($resn =~ /(DA|DG)/ && $aname =~ /N1/) || 
	     ($resn =~ /(DT|DC)/ && $aname =~ /N3/) )
	{
	    #print STDERR "N1N3: $resn.$aname\n";
	    $vstrands[$ih]{"x$dir"}[$ib][$il] = $xyz[0];
	    $vstrands[$ih]{"y$dir"}[$ib][$il] = $xyz[1];
	    $vstrands[$ih]{"z$dir"}[$ib][$il] = $xyz[2];
	}
	
	### ENM atoms
	#my $hi = which_enm_atom($resn,$aname);
	#if ($hi != -1) {
	#    $enm_atoms[$hi] = $global_aid;
	#}
	
	# H3T
	if ($resn =~ /D.3/ && $aname =~ /O3'/) {
	    ++$global_aid;
	    $atoms{$aname} = $global_aid;
	    printf "%-6s%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n",
    		"ATOM",0, " H3T", $resn, num2chr($ch), $resid, 
		$xyz[0]+rand, $xyz[1]+rand, $xyz[2]+rand,1,0, $seg;
	}
    }

    return %atoms;
}

sub num2chr {
    my $n = shift;
    $n -= 26 while ($n >= 91);
    return chr($n);
}

sub which_enm_atom {
    my $resn = shift;
    my $aname = shift;

    return 0 if ($resn =~ /DA/ && $aname eq " N6 ");
    return 0 if ($resn =~ /DT/ && $aname eq " O4 ");
    return 0 if ($resn =~ /DG/ && $aname eq " O6 ");
    return 0 if ($resn =~ /DC/ && $aname eq " N4 ");
    return 1 if ($resn =~ /DA/ && $aname eq " H2 ");
    return 1 if ($resn =~ /DT/ && $aname eq " O2 ");
    return 1 if ($resn =~ /DG/ && $aname eq " N2 ");
    return 1 if ($resn =~ /DC/ && $aname eq " O2 ");
    return 2 if ($resn =~ /(DA|DG)/ && $aname eq " N1 ");
    return 2 if ($resn =~ /(DT|DC)/ && $aname eq " N3 ");
    return 3 if ($aname eq " C1'");
    for my $i (0 .. $#baseatoms) {
	return $i+4 if ($aname eq $baseatoms[$i]);
    }
    return -1;
}

sub wc {
    my $seq = $_[0];
    die "wc(): wrong nucleotide '$seq'.\n" unless ($seq =~ /[atgcATGC]/);
    $seq =~ tr/a-z/A-Z/;
    $seq =~ tr/ATGC/TACG/;
    return $seq;
}

sub init_pdb_charmm {
   %pdb = (
    A => [
            [" P  ", [   0.288, -9.220,  -1.848 ]],
	    [" C4'", [   3.212, -6.864,  -1.355 ]],
	    [" H4'", [   4.069, -6.713,  -1.848 ]],
	    [" O4'", [   2.387, -5.664,  -1.352 ]],
	    [" C1'", [   2.281, -5.198,  -0.016 ]],
	    [" H1'", [   3.030, -4.578,   0.218 ]],
	    [" C2'", [   2.304, -6.454,   0.850 ]],
	    [" H2'", [   1.451, -6.973,   0.796 ]],
	    ["H2''", [   2.514, -6.256,   1.807 ]],
	    [" O1P", [   0.421,-10.485,  -2.605 ]],
	    [" O2P", [  -0.692, -9.226,  -0.740 ]],
	    [" O5'", [   1.721, -8.770,  -1.295 ]],
	    [" C5'", [   2.585, -7.995,  -2.146 ]],
	    [" H5'", [   2.038, -7.632,  -2.900 ]],
	    ["H5''", [   3.293, -8.601,  -2.509 ]],
	    [" N9 ", [   1.021, -4.412,   0.101 ]],
	    [" C5 ", [  -0.409, -2.735,   0.190 ]],
	    [" N7 ", [  -1.161, -3.897,   0.354 ]],
	    [" C8 ", [  -0.268, -4.856,   0.294 ]],
	    [" H8 ", [  -0.503, -5.824,   0.381 ]],
	    [" N1 ", [   0.269, -0.516,  -0.025 ]],
	    [" C2 ", [   1.584, -0.914,  -0.172 ]],
	    [" H2 ", [   2.282, -0.210,  -0.302 ]],
	    [" N3 ", [   1.969, -2.189,  -0.149 ]],
	    [" C4 ", [   0.923, -3.037,   0.035 ]],
	    [" C6 ", [  -0.820, -1.378,   0.165 ]],
	    [" N6 ", [  -1.951, -0.913,   0.284 ]],
	    [" H61", [  -2.731, -1.524,   0.420 ]],
	    [" H62", [  -2.091,  0.076,   0.245 ]],
	    [" C3'", [   3.458, -7.146,   0.127 ]],
	    [" H3'", [   3.568, -8.115,   0.350 ]],
	    [" O3'", [   4.677, -6.520,   0.518 ]]],
    T => [  
	    [" P  ", [   0.249, -9.221,  -1.866 ]],
	    [" C4'", [   3.176, -6.866,  -1.373 ]],
	    [" H4'", [   4.034, -6.717,  -1.864 ]],
	    [" O4'", [   2.348, -5.666,  -1.370 ]],
	    [" C1'", [   2.243, -5.199,  -0.034 ]],
	    [" H1'", [   2.994, -4.581,   0.200 ]],
	    [" C2'", [   2.265, -6.453,   0.832 ]],
	    [" H2'", [   1.411, -6.971,   0.778 ]],
	    ["H2''", [   2.474, -6.254,   1.789 ]],
	    [" O1P", [   0.383,-10.486,  -2.623 ]],
	    [" O2P", [  -0.730, -9.227,  -0.758 ]],
	    [" O5'", [   1.683, -8.771,  -1.313 ]],
	    [" C5'", [   2.547, -7.996,  -2.164 ]],
	    [" H5'", [   2.000, -7.632,  -2.918 ]],
	    ["H5''", [   3.254, -8.602,  -2.528 ]],
	    [" N1 ", [   0.982, -4.412,   0.083 ]],
	    [" C6 ", [  -0.215, -5.030,   0.274 ]],
	    [" H6 ", [  -0.250, -6.027,   0.336 ]],
	    [" C2 ", [   1.072, -3.026,  -0.008 ]],
	    [" O2 ", [   2.183, -2.509,  -0.181 ]],
	    [" N3 ", [  -0.065, -2.290,   0.097 ]],
	    [" H3 ", [  -0.018, -1.293,   0.034 ]],
	    [" C4 ", [  -1.249, -2.888,   0.284 ]],
	    [" O4 ", [  -2.329, -2.129,   0.380 ]],
	    [" C5 ", [  -1.361, -4.314,   0.381 ]],
	    [" C5M", [  -2.245, -4.760,   0.521 ]],
	    [" H51", [  -2.964, -4.067,   0.574 ]],
	    [" H52", [  -2.432, -5.376,  -0.244 ]],
	    [" H53", [  -2.223, -5.283,   1.373 ]],
	    [" C3'", [   3.421, -7.147,   0.109 ]],
	    [" H3'", [   3.532, -8.115,   0.332 ]],
	    [" O3'", [   4.638, -6.523,   0.500 ]]],
    G => [
	    [" P  ", [   0.288, -9.220,  -1.848 ]],
	    [" C4'", [   3.212, -6.864,  -1.355 ]],
	    [" H4'", [   4.069, -6.713,  -1.848 ]],
	    [" O4'", [   2.387, -5.664,  -1.352 ]],
	    [" C1'", [   2.281, -5.198,  -0.016 ]],
	    [" H1'", [   3.030, -4.578,   0.218 ]],
	    [" C2'", [   2.304, -6.454,   0.850 ]],
	    [" H2'", [   1.451, -6.973,   0.796 ]],
	    ["H2''", [   2.514, -6.256,   1.807 ]],
	    [" O1P", [   0.421,-10.485,  -2.605 ]],
	    [" O2P", [  -0.692, -9.226,  -0.740 ]],
	    [" O5'", [   1.721, -8.770,  -1.295 ]],
	    [" C5'", [   2.585, -7.995,  -2.146 ]],
	    [" H5'", [   2.038, -7.632,  -2.900 ]],
	    ["H5''", [   3.293, -8.601,  -2.509 ]],
	    [" N9 ", [   1.021, -4.412,   0.101 ]],
	    [" C4 ", [   0.923, -3.037,   0.035 ]],
	    [" N2 ", [   2.478,  0.061,  -0.344 ]],
	    [" H21", [   2.180,  1.015,  -0.361 ]],
	    [" H22", [   3.446, -0.163,  -0.456 ]],
	    [" N3 ", [   1.969, -2.189,  -0.149 ]],
	    [" C2 ", [   1.584, -0.914,  -0.172 ]],
	    [" N1 ", [   0.269, -0.516,  -0.025 ]],
	    [" H1 ", [   0.079,  0.465,  -0.057 ]],
	    [" C6 ", [  -0.820, -1.378,   0.165 ]],
	    [" O6 ", [  -1.951, -0.913,   0.284 ]],
	    [" C5 ", [  -0.409, -2.735,   0.190 ]],
	    [" N7 ", [  -1.161, -3.897,   0.354 ]],
	    [" C8 ", [  -0.268, -4.856,   0.294 ]],
	    [" H8 ", [  -0.503, -5.824,   0.381 ]],
	    [" C3'", [   3.458, -7.146,   0.127 ]],
	    [" H3'", [   3.568, -8.115,   0.350 ]],
	    [" O3'", [   4.677, -6.520,   0.518 ]]],
    C => [
	    [" P  ", [   0.249, -9.221,  -1.866 ]],
	    [" C4'", [   3.176, -6.866,  -1.373 ]],
	    [" H4'", [   4.034, -6.717,  -1.864 ]],
	    [" O4'", [   2.348, -5.666,  -1.370 ]],
	    [" C1'", [   2.243, -5.199,  -0.034 ]],
	    [" H1'", [   2.994, -4.581,   0.200 ]],
	    [" C2'", [   2.265, -6.453,   0.832 ]],
	    [" H2'", [   1.411, -6.971,   0.778 ]],
	    ["H2''", [   2.474, -6.254,   1.789 ]],
	    [" O1P", [   0.383,-10.486,  -2.623 ]],
	    [" O2P", [  -0.730, -9.227,  -0.758 ]],
	    [" O5'", [   1.683, -8.771,  -1.313 ]],
	    [" C5'", [   2.547, -7.996,  -2.164 ]],
	    [" H5'", [   2.000, -7.632,  -2.918 ]],
	    ["H5''", [   3.254, -8.602,  -2.528 ]],
	    [" N1 ", [   0.982, -4.412,   0.083 ]],
	    [" C6 ", [  -0.215, -5.030,   0.274 ]],
	    [" H6 ", [  -0.250, -6.027,   0.336 ]],
	    [" C5 ", [  -1.361, -4.314,   0.381 ]],
	    [" H5 ", [  -2.245, -4.760,   0.521 ]],
	    [" C2 ", [   1.072, -3.026,  -0.008 ]],
	    [" O2 ", [   2.183, -2.509,  -0.181 ]],
	    [" N3 ", [  -0.065, -2.290,   0.097 ]],
	    [" C4 ", [  -1.249, -2.888,   0.284 ]],
	    [" N4 ", [  -2.329, -2.129,   0.380 ]],
	    [" H41", [  -3.225, -2.550,   0.520 ]],
	    [" H42", [  -2.249, -1.135,   0.312 ]],
	    [" C3'", [   3.421, -7.147,   0.109 ]],
	    [" H3'", [   3.532, -8.115,   0.332 ]],
	    [" O3'", [   4.638, -6.523,   0.500 ]]],
    A3 => [
            [" P  ", [   0.288, -9.220,  -1.848 ]],
	    [" C4'", [   3.212, -6.864,  -1.355 ]],
	    [" H4'", [   4.069, -6.713,  -1.848 ]],
	    [" O4'", [   2.387, -5.664,  -1.352 ]],
	    [" C1'", [   2.281, -5.198,  -0.016 ]],
	    [" H1'", [   3.030, -4.578,   0.218 ]],
	    [" C2'", [   2.304, -6.454,   0.850 ]],
	    [" H2'", [   1.451, -6.973,   0.796 ]],
	    ["H2''", [   2.514, -6.256,   1.807 ]],
	    [" C3'", [   3.458, -7.146,   0.127 ]],
	    [" H3'", [   3.568, -8.115,   0.350 ]],
	    [" O3'", [   4.677, -6.520,   0.518 ]],
	    [" O1P", [   0.421,-10.485,  -2.605 ]],
	    [" O2P", [  -0.692, -9.226,  -0.740 ]],
	    [" O5'", [   1.721, -8.770,  -1.295 ]],
	    [" C5'", [   2.585, -7.995,  -2.146 ]],
	    [" H5'", [   2.038, -7.632,  -2.900 ]],
	    ["H5''", [   3.293, -8.601,  -2.509 ]],
	    [" N9 ", [   1.021, -4.412,   0.101 ]],
	    [" C5 ", [  -0.409, -2.735,   0.190 ]],
	    [" N7 ", [  -1.161, -3.897,   0.354 ]],
	    [" C8 ", [  -0.268, -4.856,   0.294 ]],
	    [" H8 ", [  -0.503, -5.824,   0.381 ]],
	    [" N1 ", [   0.269, -0.516,  -0.025 ]],
	    [" C2 ", [   1.584, -0.914,  -0.172 ]],
	    [" H2 ", [   2.282, -0.210,  -0.302 ]],
	    [" N3 ", [   1.969, -2.189,  -0.149 ]],
	    [" C4 ", [   0.923, -3.037,   0.035 ]],
	    [" C6 ", [  -0.820, -1.378,   0.165 ]],
	    [" N6 ", [  -1.951, -0.913,   0.284 ]],
	    [" H61", [  -2.731, -1.524,   0.420 ]],
	    [" H62", [  -2.091,  0.076,   0.245 ]]
	    ],
    T3 => [  
	    [" P  ", [   0.249, -9.221,  -1.866 ]],
	    [" C4'", [   3.176, -6.866,  -1.373 ]],
	    [" H4'", [   4.034, -6.717,  -1.864 ]],
	    [" O4'", [   2.348, -5.666,  -1.370 ]],
	    [" C1'", [   2.243, -5.199,  -0.034 ]],
	    [" H1'", [   2.994, -4.581,   0.200 ]],
	    [" C2'", [   2.265, -6.453,   0.832 ]],
	    [" H2'", [   1.411, -6.971,   0.778 ]],
	    ["H2''", [   2.474, -6.254,   1.789 ]],
	    [" C3'", [   3.421, -7.147,   0.109 ]],
	    [" H3'", [   3.532, -8.115,   0.332 ]],
	    [" O3'", [   4.638, -6.523,   0.500 ]],
	    [" O1P", [   0.383,-10.486,  -2.623 ]],
	    [" O2P", [  -0.730, -9.227,  -0.758 ]],
	    [" O5'", [   1.683, -8.771,  -1.313 ]],
	    [" C5'", [   2.547, -7.996,  -2.164 ]],
	    [" H5'", [   2.000, -7.632,  -2.918 ]],
	    ["H5''", [   3.254, -8.602,  -2.528 ]],
	    [" N1 ", [   0.982, -4.412,   0.083 ]],
	    [" C6 ", [  -0.215, -5.030,   0.274 ]],
	    [" H6 ", [  -0.250, -6.027,   0.336 ]],
	    [" C2 ", [   1.072, -3.026,  -0.008 ]],
	    [" O2 ", [   2.183, -2.509,  -0.181 ]],
	    [" N3 ", [  -0.065, -2.290,   0.097 ]],
	    [" H3 ", [  -0.018, -1.293,   0.034 ]],
	    [" C4 ", [  -1.249, -2.888,   0.284 ]],
	    [" O4 ", [  -2.329, -2.129,   0.380 ]],
	    [" C5 ", [  -1.361, -4.314,   0.381 ]],
	    [" C5M", [  -2.245, -4.760,   0.521 ]],
	    [" H51", [  -2.964, -4.067,   0.574 ]],
	    [" H52", [  -2.432, -5.376,  -0.244 ]],
	    [" H53", [  -2.223, -5.283,   1.373 ]]
	    ],
    G3 => [
	    [" P  ", [   0.288, -9.220,  -1.848 ]],
	    [" C4'", [   3.212, -6.864,  -1.355 ]],
	    [" H4'", [   4.069, -6.713,  -1.848 ]],
	    [" O4'", [   2.387, -5.664,  -1.352 ]],
	    [" C1'", [   2.281, -5.198,  -0.016 ]],
	    [" H1'", [   3.030, -4.578,   0.218 ]],
	    [" C2'", [   2.304, -6.454,   0.850 ]],
	    [" H2'", [   1.451, -6.973,   0.796 ]],
	    ["H2''", [   2.514, -6.256,   1.807 ]],
	    [" C3'", [   3.458, -7.146,   0.127 ]],
	    [" H3'", [   3.568, -8.115,   0.350 ]],
	    [" O3'", [   4.677, -6.520,   0.518 ]],
	    [" O1P", [   0.421,-10.485,  -2.605 ]],
	    [" O2P", [  -0.692, -9.226,  -0.740 ]],
	    [" O5'", [   1.721, -8.770,  -1.295 ]],
	    [" C5'", [   2.585, -7.995,  -2.146 ]],
	    [" H5'", [   2.038, -7.632,  -2.900 ]],
	    ["H5''", [   3.293, -8.601,  -2.509 ]],
	    [" N9 ", [   1.021, -4.412,   0.101 ]],
	    [" C4 ", [   0.923, -3.037,   0.035 ]],
	    [" N2 ", [   2.478,  0.061,  -0.344 ]],
	    [" H21", [   2.180,  1.015,  -0.361 ]],
	    [" H22", [   3.446, -0.163,  -0.456 ]],
	    [" N3 ", [   1.969, -2.189,  -0.149 ]],
	    [" C2 ", [   1.584, -0.914,  -0.172 ]],
	    [" N1 ", [   0.269, -0.516,  -0.025 ]],
	    [" H1 ", [   0.079,  0.465,  -0.057 ]],
	    [" C6 ", [  -0.820, -1.378,   0.165 ]],
	    [" O6 ", [  -1.951, -0.913,   0.284 ]],
	    [" C5 ", [  -0.409, -2.735,   0.190 ]],
	    [" N7 ", [  -1.161, -3.897,   0.354 ]],
	    [" C8 ", [  -0.268, -4.856,   0.294 ]],
	    [" H8 ", [  -0.503, -5.824,   0.381 ]]
	    ],
    C3 => [
	    [" P  ", [   0.249, -9.221,  -1.866 ]],
	    [" C4'", [   3.176, -6.866,  -1.373 ]],
	    [" H4'", [   4.034, -6.717,  -1.864 ]],
	    [" O4'", [   2.348, -5.666,  -1.370 ]],
	    [" C1'", [   2.243, -5.199,  -0.034 ]],
	    [" H1'", [   2.994, -4.581,   0.200 ]],
	    [" C2'", [   2.265, -6.453,   0.832 ]],
	    [" H2'", [   1.411, -6.971,   0.778 ]],
	    ["H2''", [   2.474, -6.254,   1.789 ]],
	    [" C3'", [   3.421, -7.147,   0.109 ]],
	    [" H3'", [   3.532, -8.115,   0.332 ]],
	    [" O3'", [   4.638, -6.523,   0.500 ]],
	    [" O1P", [   0.383,-10.486,  -2.623 ]],
	    [" O2P", [  -0.730, -9.227,  -0.758 ]],
	    [" O5'", [   1.683, -8.771,  -1.313 ]],
	    [" C5'", [   2.547, -7.996,  -2.164 ]],
	    [" H5'", [   2.000, -7.632,  -2.918 ]],
	    ["H5''", [   3.254, -8.602,  -2.528 ]],
	    [" N1 ", [   0.982, -4.412,   0.083 ]],
	    [" C6 ", [  -0.215, -5.030,   0.274 ]],
	    [" H6 ", [  -0.250, -6.027,   0.336 ]],
	    [" C5 ", [  -1.361, -4.314,   0.381 ]],
	    [" H5 ", [  -2.245, -4.760,   0.521 ]],
	    [" C2 ", [   1.072, -3.026,  -0.008 ]],
	    [" O2 ", [   2.183, -2.509,  -0.181 ]],
	    [" N3 ", [  -0.065, -2.290,   0.097 ]],
	    [" C4 ", [  -1.249, -2.888,   0.284 ]],
	    [" N4 ", [  -2.329, -2.129,   0.380 ]],
	    [" H41", [  -3.225, -2.550,   0.520 ]],
	    [" H42", [  -2.249, -1.135,   0.312 ]]
	    ]
    );
}

sub init_pdb_amber {
   %pdb = (
    A => [[" P  ", [   0.288, -9.220,  -1.848 ]],
	    [" O1P", [   0.421,-10.485,  -2.605 ]],
	    [" O2P", [  -0.692, -9.226,  -0.740 ]],
	    [" O5'", [   1.721, -8.770,  -1.295 ]],
	    [" C5'", [   2.585, -7.995,  -2.146 ]],
	    ["H5'1", [   2.038, -7.632,  -2.900 ]],
	    ["H5'2", [   3.293, -8.601,  -2.509 ]],
	    [" C4'", [   3.212, -6.864,  -1.355 ]],
	    [" H4'", [   4.069, -6.713,  -1.848 ]],
	    [" O4'", [   2.387, -5.664,  -1.352 ]],
	    [" C1'", [   2.281, -5.198,  -0.016 ]],
	    [" H1'", [   3.030, -4.578,   0.218 ]],
	    [" N9 ", [   1.021, -4.412,   0.101 ]],
	    [" C8 ", [  -0.268, -4.856,   0.294 ]],
	    [" H8 ", [  -0.503, -5.824,   0.381 ]],
	    [" N7 ", [  -1.161, -3.897,   0.354 ]],
	    [" C5 ", [  -0.409, -2.735,   0.190 ]],
	    [" C6 ", [  -0.820, -1.378,   0.165 ]],
	    [" N6 ", [  -1.951, -0.913,   0.284 ]],
	    [" H61", [  -2.731, -1.524,   0.420 ]],
	    [" H62", [  -2.091,  0.076,   0.245 ]],
	    [" N1 ", [   0.269, -0.516,  -0.025 ]],
	    [" C2 ", [   1.584, -0.914,  -0.172 ]],
	    [" H2 ", [   2.282, -0.210,  -0.302 ]],
	    [" N3 ", [   1.969, -2.189,  -0.149 ]],
	    [" C4 ", [   0.923, -3.037,   0.035 ]],
	    [" C3'", [   3.458, -7.146,   0.127 ]],
	    [" H3'", [   3.568, -8.115,   0.350 ]],
	    [" C2'", [   2.304, -6.454,   0.850 ]],
	    ["H2'1", [   1.451, -6.973,   0.796 ]],
	    ["H2'2", [   2.514, -6.256,   1.807 ]],
	    [" O3'", [   4.677, -6.520,   0.518 ]]],
    T => [[" P  ", [   0.249, -9.221,  -1.866 ]],
	    [" O1P", [   0.383,-10.486,  -2.623 ]],
	    [" O2P", [  -0.730, -9.227,  -0.758 ]],
	    [" O5'", [   1.683, -8.771,  -1.313 ]],
	    [" C5'", [   2.547, -7.996,  -2.164 ]],
	    ["H5'1", [   2.000, -7.632,  -2.918 ]],
	    ["H5'2", [   3.254, -8.602,  -2.528 ]],
	    [" C4'", [   3.176, -6.866,  -1.373 ]],
	    [" H4'", [   4.034, -6.717,  -1.864 ]],
	    [" O4'", [   2.348, -5.666,  -1.370 ]],
	    [" C1'", [   2.243, -5.199,  -0.034 ]],
	    [" H1'", [   2.994, -4.581,   0.200 ]],
	    [" N1 ", [   0.982, -4.412,   0.083 ]],
	    [" C6 ", [  -0.215, -5.030,   0.274 ]],
	    [" H6 ", [  -0.250, -6.027,   0.336 ]],
	    [" C5 ", [  -1.361, -4.314,   0.381 ]],
	    [" C7 ", [  -2.245, -4.760,   0.521 ]],
	    [" H71", [  -2.964, -4.067,   0.574 ]],
	    [" H72", [  -2.432, -5.376,  -0.244 ]],
	    [" H73", [  -2.223, -5.283,   1.373 ]],
	    [" C4 ", [  -1.249, -2.888,   0.284 ]],
	    [" O4 ", [  -2.329, -2.129,   0.380 ]],
	    [" N3 ", [  -0.065, -2.290,   0.097 ]],
	    [" H3 ", [  -0.018, -1.293,   0.034 ]],
	    [" C2 ", [   1.072, -3.026,  -0.008 ]],
	    [" O2 ", [   2.183, -2.509,  -0.181 ]],
	    [" C3'", [   3.421, -7.147,   0.109 ]],
	    [" H3'", [   3.532, -8.115,   0.332 ]],
	    [" C2'", [   2.265, -6.453,   0.832 ]],
	    ["H2'1", [   1.411, -6.971,   0.778 ]],
	    ["H2'2", [   2.474, -6.254,   1.789 ]],
	    [" O3'", [   4.638, -6.523,   0.500 ]]],
    G => [[" P  ", [   0.288, -9.220,  -1.848 ]],
	    [" O1P", [   0.421,-10.485,  -2.605 ]],
	    [" O2P", [  -0.692, -9.226,  -0.740 ]],
	    [" O5'", [   1.721, -8.770,  -1.295 ]],
	    [" C5'", [   2.585, -7.995,  -2.146 ]],
	    ["H5'1", [   2.038, -7.632,  -2.900 ]],
	    ["H5'2", [   3.293, -8.601,  -2.509 ]],
	    [" C4'", [   3.212, -6.864,  -1.355 ]],
	    [" H4'", [   4.069, -6.713,  -1.848 ]],
	    [" O4'", [   2.387, -5.664,  -1.352 ]],
	    [" C1'", [   2.281, -5.198,  -0.016 ]],
	    [" H1'", [   3.030, -4.578,   0.218 ]],
	    [" N9 ", [   1.021, -4.412,   0.101 ]],
	    [" C8 ", [  -0.268, -4.856,   0.294 ]],
	    [" H8 ", [  -0.503, -5.824,   0.381 ]],
	    [" N7 ", [  -1.161, -3.897,   0.354 ]],
	    [" C5 ", [  -0.409, -2.735,   0.190 ]],
	    [" C6 ", [  -0.820, -1.378,   0.165 ]],
	    [" O6 ", [  -1.951, -0.913,   0.284 ]],
	    [" N1 ", [   0.269, -0.516,  -0.025 ]],
	    [" H1 ", [   0.079,  0.465,  -0.057 ]],
	    [" C2 ", [   1.584, -0.914,  -0.172 ]],
	    [" N2 ", [   2.478,  0.061,  -0.344 ]],
	    [" H21", [   2.180,  1.015,  -0.361 ]],
	    [" H22", [   3.446, -0.163,  -0.456 ]],
	    [" N3 ", [   1.969, -2.189,  -0.149 ]],
	    [" C4 ", [   0.923, -3.037,   0.035 ]],
	    [" C3'", [   3.458, -7.146,   0.127 ]],
	    [" H3'", [   3.568, -8.115,   0.350 ]],
	    [" C2'", [   2.304, -6.454,   0.850 ]],
	    ["H2'1", [   1.451, -6.973,   0.796 ]],
	    ["H2'2", [   2.514, -6.256,   1.807 ]],
	    [" O3'", [   4.677, -6.520,   0.518 ]]],
    C => [[" P  ", [   0.249, -9.221,  -1.866 ]],
	    [" O1P", [   0.383,-10.486,  -2.623 ]],
	    [" O2P", [  -0.730, -9.227,  -0.758 ]],
	    [" O5'", [   1.683, -8.771,  -1.313 ]],
	    [" C5'", [   2.547, -7.996,  -2.164 ]],
	    ["H5'1", [   2.000, -7.632,  -2.918 ]],
	    ["H5'2", [   3.254, -8.602,  -2.528 ]],
	    [" C4'", [   3.176, -6.866,  -1.373 ]],
	    [" H4'", [   4.034, -6.717,  -1.864 ]],
	    [" O4'", [   2.348, -5.666,  -1.370 ]],
	    [" C1'", [   2.243, -5.199,  -0.034 ]],
	    [" H1'", [   2.994, -4.581,   0.200 ]],
	    [" N1 ", [   0.982, -4.412,   0.083 ]],
	    [" C6 ", [  -0.215, -5.030,   0.274 ]],
	    [" H6 ", [  -0.250, -6.027,   0.336 ]],
	    [" C5 ", [  -1.361, -4.314,   0.381 ]],
	    [" H5 ", [  -2.245, -4.760,   0.521 ]],
	    [" C4 ", [  -1.249, -2.888,   0.284 ]],
	    [" N4 ", [  -2.329, -2.129,   0.380 ]],
	    [" H41", [  -3.225, -2.550,   0.520 ]],
	    [" H42", [  -2.249, -1.135,   0.312 ]],
	    [" N3 ", [  -0.065, -2.290,   0.097 ]],
	    [" C2 ", [   1.072, -3.026,  -0.008 ]],
	    [" O2 ", [   2.183, -2.509,  -0.181 ]],
	    [" C3'", [   3.421, -7.147,   0.109 ]],
	    [" H3'", [   3.532, -8.115,   0.332 ]],
	    [" C2'", [   2.265, -6.453,   0.832 ]],
	    ["H2'1", [   1.411, -6.971,   0.778 ]],
	    ["H2'2", [   2.474, -6.254,   1.789 ]],
	    [" O3'", [   4.638, -6.523,   0.500 ]]]
    );
}

sub flip_x {
    return ($_[0], -$_[1], -$_[2]);
}

sub rotate_z {
    my $x = shift;
    my $y = shift;
    my $z = shift;
    my $th = shift;
    $th = $th/180*3.14159265359;

    return ($x*cos($th)-$y*sin($th), $x*sin($th)+$y*cos($th), $z);
}

sub vec_sub {
    return ($_[0]-$_[3],$_[1]-$_[4],$_[2]-$_[5]);
}

sub vec_add {
    return ($_[0]+$_[3],$_[1]+$_[4],$_[2]+$_[5]);
}

sub vec_norm {
    return sqrt($_[0]*$_[0]+$_[1]*$_[1]+$_[2]*$_[2]);
}

sub vec_unit {
    my $len = vec_norm(@_);

    return ($_[0]/$len,$_[1]/$len,$_[2]/$len);
}

sub vec_inner {
    return ($_[0]*$_[3]+$_[1]*$_[4]+$_[2]*$_[5]);
}

sub vec_cross {
    return ($_[1]*$_[5]-$_[2]*$_[4],
            $_[2]*$_[3]-$_[0]*$_[5],
            $_[0]*$_[4]-$_[1]*$_[3]);
}

sub na1to3 {
    if ($_[0] eq 'A') {
	return "DA";
    }
    elsif ($_[0] eq 'T') {
	return "DT";
    }
    elsif ($_[0] eq 'G') {
	return "DG";
    }
    elsif ($_[0] eq 'C') {
	return "DC";
    }
}

sub init_enm {
    $enm_stack{DADA}{' C6 '}{' C6 '} = 3.51952582033432;
    $enm_stack{DADA}{' C6 '}{' N1 '} = 3.603382438765;
    $enm_stack{DADA}{' C6 '}{' N9 '} = 5.58112927640993;
    $enm_stack{DADA}{' C6 '}{' C5 '} = 4.12637383182862;
    $enm_stack{DADA}{' C6 '}{' N3 '} = 4.80445668104105;
    $enm_stack{DADA}{' C6 '}{' N7 '} = 4.83621049996792;
    $enm_stack{DADA}{' C6 '}{' C4 '} = 4.67675453706948;
    $enm_stack{DADA}{' C6 '}{' C8 '} = 5.60293931789378;
    $enm_stack{DADA}{' C6 '}{' C2 '} = 4.22934616696245;
    $enm_stack{DADA}{' N1 '}{' C6 '} = 3.72475784447795;
    $enm_stack{DADA}{' N1 '}{' N1 '} = 3.3920790085138;
    $enm_stack{DADA}{' N1 '}{' N9 '} = 5.36425437502734;
    $enm_stack{DADA}{' N1 '}{' C5 '} = 4.24052048692139;
    $enm_stack{DADA}{' N1 '}{' N3 '} = 4.18331029210122;
    $enm_stack{DADA}{' N1 '}{' N7 '} = 5.17133106656304;
    $enm_stack{DADA}{' N1 '}{' C4 '} = 4.392304178902;
    $enm_stack{DADA}{' N1 '}{' C8 '} = 5.71328714489303;
    $enm_stack{DADA}{' N1 '}{' C2 '} = 3.62672110866;
    $enm_stack{DADA}{' N9 '}{' C6 '} = 4.52605357900235;
    $enm_stack{DADA}{' N9 '}{' N1 '} = 5.29648015194997;
    $enm_stack{DADA}{' N9 '}{' N9 '} = 4.38567976487112;
    $enm_stack{DADA}{' N9 '}{' C5 '} = 3.98200426920916;
    $enm_stack{DADA}{' N9 '}{' N3 '} = 5.23927504527105;
    $enm_stack{DADA}{' N9 '}{' N7 '} = 3.69051364988669;
    $enm_stack{DADA}{' N9 '}{' C4 '} = 4.41336810157503;
    $enm_stack{DADA}{' N9 '}{' C8 '} = 3.93670141108009;
    $enm_stack{DADA}{' N9 '}{' C2 '} = 5.53976948978926;
    $enm_stack{DADA}{' C5 '}{' C6 '} = 3.59375959129155;
    $enm_stack{DADA}{' C5 '}{' N1 '} = 4.110139778645;
    $enm_stack{DADA}{' C5 '}{' N9 '} = 5.05601681168091;
    $enm_stack{DADA}{' C5 '}{' C5 '} = 3.78928634441896;
    $enm_stack{DADA}{' C5 '}{' N3 '} = 4.93103498263803;
    $enm_stack{DADA}{' C5 '}{' N7 '} = 4.11286445193614;
    $enm_stack{DADA}{' C5 '}{' C4 '} = 4.44836205810633;
    $enm_stack{DADA}{' C5 '}{' C8 '} = 4.8286904021691;
    $enm_stack{DADA}{' C5 '}{' C2 '} = 4.68213850713539;
    $enm_stack{DADA}{' N3 '}{' C6 '} = 4.11427988838873;
    $enm_stack{DADA}{' N3 '}{' N1 '} = 4.2462616499693;
    $enm_stack{DADA}{' N3 '}{' N9 '} = 4.00042810209107;
    $enm_stack{DADA}{' N3 '}{' C5 '} = 3.78168057879033;
    $enm_stack{DADA}{' N3 '}{' N3 '} = 3.82054040680111;
    $enm_stack{DADA}{' N3 '}{' N7 '} = 4.26424541976655;
    $enm_stack{DADA}{' N3 '}{' C4 '} = 3.61687904138361;
    $enm_stack{DADA}{' N3 '}{' C8 '} = 4.33603782271327;
    $enm_stack{DADA}{' N3 '}{' C2 '} = 4.06822590326545;
    $enm_stack{DADA}{' N7 '}{' C6 '} = 4.16490396047736;
    $enm_stack{DADA}{' N7 '}{' N1 '} = 4.98258968810397;
    $enm_stack{DADA}{' N7 '}{' N9 '} = 5.62022499549618;
    $enm_stack{DADA}{' N7 '}{' C5 '} = 4.28311218624962;
    $enm_stack{DADA}{' N7 '}{' N3 '} = 5.91775675404118;
    $enm_stack{DADA}{' N7 '}{' N7 '} = 4.21512336711513;
    $enm_stack{DADA}{' N7 '}{' C4 '} = 5.19388794642318;
    $enm_stack{DADA}{' N7 '}{' C8 '} = 5.05215498574618;
    $enm_stack{DADA}{' N7 '}{' C2 '} = 5.71742022244298;
    $enm_stack{DADA}{' C4 '}{' C6 '} = 3.86482315248706;
    $enm_stack{DADA}{' C4 '}{' N1 '} = 4.35089979199705;
    $enm_stack{DADA}{' C4 '}{' N9 '} = 4.25125581446236;
    $enm_stack{DADA}{' C4 '}{' C5 '} = 3.59672086767934;
    $enm_stack{DADA}{' C4 '}{' N3 '} = 4.439019711603;
    $enm_stack{DADA}{' C4 '}{' N7 '} = 3.8093600512422;
    $enm_stack{DADA}{' C4 '}{' C4 '} = 3.90640051710011;
    $enm_stack{DADA}{' C4 '}{' C8 '} = 4.15856237178186;
    $enm_stack{DADA}{' C4 '}{' C2 '} = 4.55250052169135;
    $enm_stack{DADA}{' C8 '}{' C6 '} = 4.63549371696263;
    $enm_stack{DADA}{' C8 '}{' N1 '} = 5.56668833329117;
    $enm_stack{DADA}{' C8 '}{' N9 '} = 5.22397530622035;
    $enm_stack{DADA}{' C8 '}{' C5 '} = 4.34120547774463;
    $enm_stack{DADA}{' C8 '}{' N3 '} = 6.01952572882615;
    $enm_stack{DADA}{' C8 '}{' N7 '} = 3.93796470273668;
    $enm_stack{DADA}{' C8 '}{' C4 '} = 5.12430083035725;
    $enm_stack{DADA}{' C8 '}{' C8 '} = 4.52599878479878;
    $enm_stack{DADA}{' C8 '}{' C2 '} = 6.10689487382909;
    $enm_stack{DADA}{' C2 '}{' C6 '} = 3.96438923416962;
    $enm_stack{DADA}{' C2 '}{' N1 '} = 3.71303366534698;
    $enm_stack{DADA}{' C2 '}{' N9 '} = 4.61205225469096;
    $enm_stack{DADA}{' C2 '}{' C5 '} = 4.03943968886775;
    $enm_stack{DADA}{' C2 '}{' N3 '} = 3.66725701308212;
    $enm_stack{DADA}{' C2 '}{' N7 '} = 4.86976159580733;
    $enm_stack{DADA}{' C2 '}{' C4 '} = 3.85995686504396;
    $enm_stack{DADA}{' C2 '}{' C8 '} = 5.11358670602152;
    $enm_stack{DADA}{' C2 '}{' C2 '} = 3.5400538131503;

    $enm_stack{DADT}{' C6 '}{' N1 '} = 5.57666620840803;
    $enm_stack{DADT}{' C6 '}{' C6 '} = 5.72380249484554;
    $enm_stack{DADT}{' C6 '}{' C5 '} = 5.15718392148273;
    $enm_stack{DADT}{' C6 '}{' N3 '} = 3.9918678585344;
    $enm_stack{DADT}{' C6 '}{' C4 '} = 4.1769609765953;
    $enm_stack{DADT}{' C6 '}{' C2 '} = 4.76306634427865;
    $enm_stack{DADT}{' N1 '}{' N1 '} = 5.35890492544886;
    $enm_stack{DADT}{' N1 '}{' C6 '} = 5.82186799919064;
    $enm_stack{DADT}{' N1 '}{' C5 '} = 5.53356431244817;
    $enm_stack{DADT}{' N1 '}{' N3 '} = 3.97086552277964;
    $enm_stack{DADT}{' N1 '}{' C4 '} = 4.56775426221683;
    $enm_stack{DADT}{' N1 '}{' C2 '} = 4.41998065606627;
    $enm_stack{DADT}{' N9 '}{' N1 '} = 4.38154185190556;
    $enm_stack{DADT}{' N9 '}{' C6 '} = 3.98446407939637;
    $enm_stack{DADT}{' N9 '}{' C5 '} = 3.70303510650385;
    $enm_stack{DADT}{' N9 '}{' N3 '} = 4.20903266796541;
    $enm_stack{DADT}{' N9 '}{' C4 '} = 3.83116405286957;
    $enm_stack{DADT}{' N9 '}{' C2 '} = 4.52032498831666;
    $enm_stack{DADT}{' C5 '}{' N1 '} = 5.05201326205702;
    $enm_stack{DADT}{' C5 '}{' C6 '} = 4.93547170997869;
    $enm_stack{DADT}{' C5 '}{' C5 '} = 4.33579635130618;
    $enm_stack{DADT}{' C5 '}{' N3 '} = 3.82931991350945;
    $enm_stack{DADT}{' C5 '}{' C4 '} = 3.67226537712078;
    $enm_stack{DADT}{' C5 '}{' C2 '} = 4.56342294774438;
    $enm_stack{DADT}{' N3 '}{' N1 '} = 3.99433411221445;
    $enm_stack{DADT}{' N3 '}{' C6 '} = 4.39386595152834;
    $enm_stack{DADT}{' N3 '}{' C5 '} = 4.49439662246224;
    $enm_stack{DADT}{' N3 '}{' N3 '} = 3.69098848548732;
    $enm_stack{DADT}{' N3 '}{' C4 '} = 4.15122271144298;
    $enm_stack{DADT}{' N3 '}{' C2 '} = 3.6274406955869;
    $enm_stack{DADT}{' N7 '}{' N1 '} = 5.61722751898123;
    $enm_stack{DADT}{' N7 '}{' C6 '} = 5.15531085386711;
    $enm_stack{DADT}{' N7 '}{' C5 '} = 4.3251451998748;
    $enm_stack{DADT}{' N7 '}{' N3 '} = 4.47301285488875;
    $enm_stack{DADT}{' N7 '}{' C4 '} = 3.87235819624167;
    $enm_stack{DADT}{' N7 '}{' C2 '} = 5.34439500785636;
    $enm_stack{DADT}{' C4 '}{' N1 '} = 4.24641707796114;
    $enm_stack{DADT}{' C4 '}{' C6 '} = 4.23354567708912;
    $enm_stack{DADT}{' C4 '}{' C5 '} = 3.98569642597125;
    $enm_stack{DADT}{' C4 '}{' N3 '} = 3.65351953600908;
    $enm_stack{DADT}{' C4 '}{' C4 '} = 3.66353995474323;
    $enm_stack{DADT}{' C4 '}{' C2 '} = 3.98821601721873;
    $enm_stack{DADT}{' C8 '}{' N1 '} = 5.22089063283268;
    $enm_stack{DADT}{' C8 '}{' C6 '} = 4.59697922553496;
    $enm_stack{DADT}{' C8 '}{' C5 '} = 3.92437065018074;
    $enm_stack{DADT}{' C8 '}{' N3 '} = 4.62511556612373;
    $enm_stack{DADT}{' C8 '}{' C4 '} = 3.92556569171884;
    $enm_stack{DADT}{' C8 '}{' C2 '} = 5.27026318128421;
    $enm_stack{DADT}{' C2 '}{' N1 '} = 4.60533190986274;
    $enm_stack{DADT}{' C2 '}{' C6 '} = 5.19434875610023;
    $enm_stack{DADT}{' C2 '}{' C5 '} = 5.19251827151335;
    $enm_stack{DADT}{' C2 '}{' N3 '} = 3.79119308397765;
    $enm_stack{DADT}{' C2 '}{' C4 '} = 4.50114329920744;
    $enm_stack{DADT}{' C2 '}{' C2 '} = 3.8491200293054;

    $enm_stack{DTDA}{' N1 '}{' C6 '} = 4.51452777153935;
    $enm_stack{DTDA}{' N1 '}{' N1 '} = 5.28418432683797;
    $enm_stack{DTDA}{' N1 '}{' N9 '} = 4.38148525046017;
    $enm_stack{DTDA}{' N1 '}{' C5 '} = 3.97329598696095;
    $enm_stack{DTDA}{' N1 '}{' N3 '} = 5.23097294965287;
    $enm_stack{DTDA}{' N1 '}{' N7 '} = 3.68499796472128;
    $enm_stack{DTDA}{' N1 '}{' C4 '} = 4.40601793005884;
    $enm_stack{DTDA}{' N1 '}{' C8 '} = 3.93486505486529;
    $enm_stack{DTDA}{' N1 '}{' C2 '} = 5.52849717373537;
    $enm_stack{DTDA}{' C6 '}{' C6 '} = 4.72566725447317;
    $enm_stack{DTDA}{' C6 '}{' N1 '} = 5.67215003327662;
    $enm_stack{DTDA}{' C6 '}{' N9 '} = 5.23459081495393;
    $enm_stack{DTDA}{' C6 '}{' C5 '} = 4.39665611573159;
    $enm_stack{DTDA}{' C6 '}{' N3 '} = 6.08877885950869;
    $enm_stack{DTDA}{' C6 '}{' N7 '} = 3.94982303401051;
    $enm_stack{DTDA}{' C6 '}{' C4 '} = 5.17346054396861;
    $enm_stack{DTDA}{' C6 '}{' C8 '} = 4.51363057859192;
    $enm_stack{DTDA}{' C6 '}{' C2 '} = 6.20051828478878;
    $enm_stack{DTDA}{' C5 '}{' C6 '} = 4.44009932321339;
    $enm_stack{DTDA}{' C5 '}{' N1 '} = 5.33177240699563;
    $enm_stack{DTDA}{' C5 '}{' N9 '} = 5.84047172752339;
    $enm_stack{DTDA}{' C5 '}{' C5 '} = 4.51834992004825;
    $enm_stack{DTDA}{' C5 '}{' N3 '} = 6.27017742970644;
    $enm_stack{DTDA}{' C5 '}{' N7 '} = 4.32444805726696;
    $enm_stack{DTDA}{' C5 '}{' C4 '} = 5.48095785424409;
    $enm_stack{DTDA}{' C5 '}{' C8 '} = 5.17499893719796;
    $enm_stack{DTDA}{' C5 '}{' C2 '} = 6.09438676160284;
    $enm_stack{DTDA}{' N3 '}{' C6 '} = 3.50553191398966;
    $enm_stack{DTDA}{' N3 '}{' N1 '} = 3.8747419268901;
    $enm_stack{DTDA}{' N3 '}{' N9 '} = 4.86974845346246;
    $enm_stack{DTDA}{' N3 '}{' C5 '} = 3.69563336926162;
    $enm_stack{DTDA}{' N3 '}{' N3 '} = 4.58152038083429;
    $enm_stack{DTDA}{' N3 '}{' N7 '} = 4.15542765548866;
    $enm_stack{DTDA}{' N3 '}{' C4 '} = 4.21021365253594;
    $enm_stack{DTDA}{' N3 '}{' C8 '} = 4.78986074954168;
    $enm_stack{DTDA}{' N3 '}{' C2 '} = 4.33647391321567;
    $enm_stack{DTDA}{' C4 '}{' C6 '} = 3.78319349227607;
    $enm_stack{DTDA}{' C4 '}{' N1 '} = 4.37973172694401;
    $enm_stack{DTDA}{' C4 '}{' N9 '} = 5.69210567365013;
    $enm_stack{DTDA}{' C4 '}{' C5 '} = 4.17937926012943;
    $enm_stack{DTDA}{' C4 '}{' N3 '} = 5.55256355209015;
    $enm_stack{DTDA}{' C4 '}{' N7 '} = 4.44332330581514;
    $enm_stack{DTDA}{' C4 '}{' C4 '} = 5.03576349722661;
    $enm_stack{DTDA}{' C4 '}{' C8 '} = 5.32998921199659;
    $enm_stack{DTDA}{' C4 '}{' C2 '} = 5.16168858030006;
    $enm_stack{DTDA}{' C2 '}{' C6 '} = 3.92529056249343;
    $enm_stack{DTDA}{' C2 '}{' N1 '} = 4.38874754343423;
    $enm_stack{DTDA}{' C2 '}{' N9 '} = 4.15484536415015;
    $enm_stack{DTDA}{' C2 '}{' C5 '} = 3.60366979619388;
    $enm_stack{DTDA}{' C2 '}{' N3 '} = 4.36812328122731;
    $enm_stack{DTDA}{' C2 '}{' N7 '} = 3.81527154996863;
    $enm_stack{DTDA}{' C2 '}{' C4 '} = 3.84498127433672;
    $enm_stack{DTDA}{' C2 '}{' C8 '} = 4.10256188253145;
    $enm_stack{DTDA}{' C2 '}{' C2 '} = 4.53091668429248;

    $enm_stack{DTDT}{' N1 '}{' N1 '} = 4.38746612522536;
    $enm_stack{DTDT}{' N1 '}{' C6 '} = 3.99186197156164;
    $enm_stack{DTDT}{' N1 '}{' C5 '} = 3.70869478388287;
    $enm_stack{DTDT}{' N1 '}{' N3 '} = 4.21339934020026;
    $enm_stack{DTDT}{' N1 '}{' C4 '} = 3.83677195048129;
    $enm_stack{DTDT}{' N1 '}{' C2 '} = 4.52445013233653;
    $enm_stack{DTDT}{' C6 '}{' N1 '} = 5.24036544527192;
    $enm_stack{DTDT}{' C6 '}{' C6 '} = 4.58920679420747;
    $enm_stack{DTDT}{' C6 '}{' C5 '} = 3.92769728976152;
    $enm_stack{DTDT}{' C6 '}{' N3 '} = 4.70399755527147;
    $enm_stack{DTDT}{' C6 '}{' C4 '} = 3.98587782050579;
    $enm_stack{DTDT}{' C6 '}{' C2 '} = 5.32865592809294;
    $enm_stack{DTDT}{' C5 '}{' N1 '} = 5.84305476613047;
    $enm_stack{DTDT}{' C5 '}{' C6 '} = 5.27884892755987;
    $enm_stack{DTDT}{' C5 '}{' C5 '} = 4.39600398088992;
    $enm_stack{DTDT}{' C5 '}{' N3 '} = 4.75864728678224;
    $enm_stack{DTDT}{' C5 '}{' C4 '} = 4.04096535496161;
    $enm_stack{DTDT}{' C5 '}{' C2 '} = 5.64569508563826;
    $enm_stack{DTDT}{' N3 '}{' N1 '} = 4.86842787766236;
    $enm_stack{DTDT}{' N3 '}{' C6 '} = 4.89640122947456;
    $enm_stack{DTDT}{' N3 '}{' C5 '} = 4.41387788684735;
    $enm_stack{DTDT}{' N3 '}{' N3 '} = 3.68071419700036;
    $enm_stack{DTDT}{' N3 '}{' C4 '} = 3.7175074983112;
    $enm_stack{DTDT}{' N3 '}{' C2 '} = 4.31139235050581;
    $enm_stack{DTDT}{' C4 '}{' N1 '} = 5.69132796454395;
    $enm_stack{DTDT}{' C4 '}{' C6 '} = 5.44907138877809;
    $enm_stack{DTDT}{' C4 '}{' C5 '} = 4.65535326264291;
    $enm_stack{DTDT}{' C4 '}{' N3 '} = 4.25726872536841;
    $enm_stack{DTDT}{' C4 '}{' C4 '} = 3.91059458394756;
    $enm_stack{DTDT}{' C4 '}{' C2 '} = 5.17574864150105;
    $enm_stack{DTDT}{' C2 '}{' N1 '} = 4.15699963916284;
    $enm_stack{DTDT}{' C2 '}{' C6 '} = 4.17794123941446;
    $enm_stack{DTDT}{' C2 '}{' C5 '} = 3.99521363634036;
    $enm_stack{DTDT}{' C2 '}{' N3 '} = 3.66746165624128;
    $enm_stack{DTDT}{' C2 '}{' C4 '} = 3.71881782291093;
    $enm_stack{DTDT}{' C2 '}{' C2 '} = 3.92873440181441;

    $enm_bp{DADT}{' C6 '}{' N1 '} = 6.03330705334976;
    $enm_bp{DADT}{' C6 '}{' C6 '} = 6.36104731942783;
    $enm_bp{DADT}{' C6 '}{' C5 '} = 5.6684242960456;
    $enm_bp{DADT}{' C6 '}{' N3 '} = 3.76655080942764;
    $enm_bp{DADT}{' C6 '}{' C4 '} = 4.22601396116956;
    $enm_bp{DADT}{' C6 '}{' C2 '} = 4.80778805689269;
    $enm_bp{DADT}{' N1 '}{' N1 '} = 4.99883216361582;
    $enm_bp{DADT}{' N1 '}{' C6 '} = 5.54201497652253;
    $enm_bp{DADT}{' N1 '}{' C5 '} = 5.08979616487733;
    $enm_bp{DADT}{' N1 '}{' N3 '} = 2.87864794651934;
    $enm_bp{DADT}{' N1 '}{' C4 '} = 3.71500134589478;
    $enm_bp{DADT}{' N1 '}{' C2 '} = 3.68892301356371;
    $enm_bp{DADT}{' N9 '}{' N1 '} = 8.94199066203941;
    $enm_bp{DADT}{' N9 '}{' C6 '} = 9.556211383179;
    $enm_bp{DADT}{' N9 '}{' C5 '} = 9.04706300409144;
    $enm_bp{DADT}{' N9 '}{' N3 '} = 6.89347416039257;
    $enm_bp{DADT}{' N9 '}{' C4 '} = 7.62594394157209;
    $enm_bp{DADT}{' N9 '}{' C2 '} = 7.58560590856129;
    $enm_bp{DADT}{' C5 '}{' N1 '} = 7.3139361495709;
    $enm_bp{DADT}{' C5 '}{' C6 '} = 7.72733078106535;
    $enm_bp{DADT}{' C5 '}{' C5 '} = 7.0671230355782;
    $enm_bp{DADT}{' C5 '}{' N3 '} = 5.0972936937163;
    $enm_bp{DADT}{' C5 '}{' C4 '} = 5.62465412269946;
    $enm_bp{DADT}{' C5 '}{' C2 '} = 6.02934440880598;
    $enm_bp{DADT}{' N3 '}{' N1 '} = 6.8364236995669;
    $enm_bp{DADT}{' N3 '}{' C6 '} = 7.63909196174519;
    $enm_bp{DADT}{' N3 '}{' C5 '} = 7.37695296175867;
    $enm_bp{DADT}{' N3 '}{' N3 '} = 5.04912982601953;
    $enm_bp{DADT}{' N3 '}{' C4 '} = 6.06223465398693;
    $enm_bp{DADT}{' N3 '}{' C2 '} = 5.46554123577894;
    $enm_bp{DADT}{' N7 '}{' N1 '} = 8.58416233537088;
    $enm_bp{DADT}{' N7 '}{' C6 '} = 8.89297981556238;
    $enm_bp{DADT}{' N7 '}{' C5 '} = 8.11695835889282;
    $enm_bp{DADT}{' N7 '}{' N3 '} = 6.32028179751504;
    $enm_bp{DADT}{' N7 '}{' C4 '} = 6.67658977323004;
    $enm_bp{DADT}{' N7 '}{' C2 '} = 7.33634929648255;
    $enm_bp{DADT}{' C4 '}{' N1 '} = 7.56750698711273;
    $enm_bp{DADT}{' C4 '}{' C6 '} = 8.19452323201295;
    $enm_bp{DADT}{' C4 '}{' C5 '} = 7.72696305154878;
    $enm_bp{DADT}{' C4 '}{' N3 '} = 5.53383348141232;
    $enm_bp{DADT}{' C4 '}{' C4 '} = 6.32211309610956;
    $enm_bp{DADT}{' C4 '}{' C2 '} = 6.21228412099769;
    $enm_bp{DADT}{' C8 '}{' N1 '} = 9.41119269805905;
    $enm_bp{DADT}{' C8 '}{' C6 '} = 9.85392962223701;
    $enm_bp{DADT}{' C8 '}{' C5 '} = 9.17697542766679;
    $enm_bp{DADT}{' C8 '}{' N3 '} = 7.21887283445276;
    $enm_bp{DADT}{' C8 '}{' C4 '} = 7.73207475390661;
    $enm_bp{DADT}{' C8 '}{' C2 '} = 8.10160095289814;
    $enm_bp{DADT}{' C2 '}{' N1 '} = 5.52138605786627;
    $enm_bp{DADT}{' C2 '}{' C6 '} = 6.30772994031926;
    $enm_bp{DADT}{' C2 '}{' C5 '} = 6.0741065186577;
    $enm_bp{DADT}{' C2 '}{' N3 '} = 3.72849420007595;
    $enm_bp{DADT}{' C2 '}{' C4 '} = 4.79525181820517;
    $enm_bp{DADT}{' C2 '}{' C2 '} = 4.15190245550157;

    $enm_bp{DTDA}{' N1 '}{' C6 '} = 6.03359312516182;
    $enm_bp{DTDA}{' N1 '}{' N1 '} = 4.99905891143523;
    $enm_bp{DTDA}{' N1 '}{' N9 '} = 8.94168602669541;
    $enm_bp{DTDA}{' N1 '}{' C5 '} = 7.31564255551076;
    $enm_bp{DTDA}{' N1 '}{' N3 '} = 6.83715796219453;
    $enm_bp{DTDA}{' N1 '}{' N7 '} = 8.58487676090927;
    $enm_bp{DTDA}{' N1 '}{' C4 '} = 7.56874163913659;
    $enm_bp{DTDA}{' N1 '}{' C8 '} = 9.41223687547227;
    $enm_bp{DTDA}{' N1 '}{' C2 '} = 5.52226927992469;
    $enm_bp{DTDA}{' C6 '}{' C6 '} = 6.36070507098073;
    $enm_bp{DTDA}{' C6 '}{' N1 '} = 5.54173447577561;
    $enm_bp{DTDA}{' C6 '}{' N9 '} = 9.55532260051957;
    $enm_bp{DTDA}{' C6 '}{' C5 '} = 7.72834115447811;
    $enm_bp{DTDA}{' C6 '}{' N3 '} = 7.63931423623875;
    $enm_bp{DTDA}{' C6 '}{' N7 '} = 8.8930742715891;
    $enm_bp{DTDA}{' C6 '}{' C4 '} = 8.19505527742187;
    $enm_bp{DTDA}{' C6 '}{' C8 '} = 9.85438044729348;
    $enm_bp{DTDA}{' C6 '}{' C2 '} = 6.3081950667366;
    $enm_bp{DTDA}{' C5 '}{' C6 '} = 5.66851400280532;
    $enm_bp{DTDA}{' C5 '}{' N1 '} = 5.08979950489211;
    $enm_bp{DTDA}{' C5 '}{' N9 '} = 9.04645378034951;
    $enm_bp{DTDA}{' C5 '}{' C5 '} = 7.06843440940071;
    $enm_bp{DTDA}{' C5 '}{' N3 '} = 7.37722373254329;
    $enm_bp{DTDA}{' C5 '}{' N7 '} = 8.11757365226827;
    $enm_bp{DTDA}{' C5 '}{' C4 '} = 7.72757562240577;
    $enm_bp{DTDA}{' C5 '}{' C8 '} = 9.17786173354121;
    $enm_bp{DTDA}{' C5 '}{' C2 '} = 6.07467225453357;
    $enm_bp{DTDA}{' N3 '}{' C6 '} = 3.76598313857086;
    $enm_bp{DTDA}{' N3 '}{' N1 '} = 2.87809954657583;
    $enm_bp{DTDA}{' N3 '}{' N9 '} = 6.89231934837613;
    $enm_bp{DTDA}{' N3 '}{' C5 '} = 5.09809307486633;
    $enm_bp{DTDA}{' N3 '}{' N3 '} = 5.04909318590972;
    $enm_bp{DTDA}{' N3 '}{' N7 '} = 6.32013654915778;
    $enm_bp{DTDA}{' N3 '}{' C4 '} = 5.53407869477838;
    $enm_bp{DTDA}{' N3 '}{' C8 '} = 7.21906849946723;
    $enm_bp{DTDA}{' N3 '}{' C2 '} = 3.72873611294766;
    $enm_bp{DTDA}{' C4 '}{' C6 '} = 4.22586925969084;
    $enm_bp{DTDA}{' C4 '}{' N1 '} = 3.71502072672549;
    $enm_bp{DTDA}{' C4 '}{' N9 '} = 7.62523848807367;
    $enm_bp{DTDA}{' C4 '}{' C5 '} = 5.62573559634649;
    $enm_bp{DTDA}{' C4 '}{' N3 '} = 6.06261247318349;
    $enm_bp{DTDA}{' C4 '}{' N7 '} = 6.67687127328362;
    $enm_bp{DTDA}{' C4 '}{' C4 '} = 6.32261567707543;
    $enm_bp{DTDA}{' C4 '}{' C8 '} = 7.73271634032957;
    $enm_bp{DTDA}{' C4 '}{' C2 '} = 4.79598780231977;
    $enm_bp{DTDA}{' C2 '}{' C6 '} = 4.80741448181868;
    $enm_bp{DTDA}{' C2 '}{' N1 '} = 3.68830448851501;
    $enm_bp{DTDA}{' C2 '}{' N9 '} = 7.5843625308921;
    $enm_bp{DTDA}{' C2 '}{' C5 '} = 6.03029816841589;
    $enm_bp{DTDA}{' C2 '}{' N3 '} = 5.46523366380615;
    $enm_bp{DTDA}{' C2 '}{' N7 '} = 7.33633607463562;
    $enm_bp{DTDA}{' C2 '}{' C4 '} = 6.21260348002349;
    $enm_bp{DTDA}{' C2 '}{' C8 '} = 8.10180603322494;
    $enm_bp{DTDA}{' C2 '}{' C2 '} = 4.1517439709115;
}

## 
# trace 5 ib's straight down the helix (no crossing)
# Return array of [ih,ib,il]
sub trace_helix_upstream {
    my $ih = shift;
    my $ib0 = shift;

    my @ret;
    for my $ib ($ib0 .. ($ib0+4)) {
	next if $vstrands[$ih]->{skip}[$ib] == -1;
	
	for my $il (0 .. $vstrands[$ih]->{loop}[$ib]) {
	    push @ret, [$ih,$ib,$il];
	}
    }
    return @ret;
}

sub trace_helix_downstream {
    my $ih = shift;
    my $ib0 = shift;

    my @ret;
    for my $ib (($ib0-4) .. $ib0) {
	next if $vstrands[$ih]->{skip}[$ib] == -1;
	
	for my $il (0 .. $vstrands[$ih]->{loop}[$ib]) {
	    push @ret, [$ih,$ib,$il];
	}
    }
    return reverse(@ret);
}

sub next_nonskip_ih_ib {
    my $stap_or_scaf = shift;
    my $ih = shift;
    my $ib = shift;
    while ($vstrands[$ih]{"skip"}[$ib] == -1) {
	$ih = $vstrands[$ih]{$stap_or_scaf}[$ib][2]; 
	$ib = $vstrands[$ih]{$stap_or_scaf}[$ib][3]; 
	#last if $ih == -1;
    }
    return ($ih,$ib);
}

sub get_direction_scaf {
    my $ih = shift;
    my $sum = $vstrands[$ih]{row} + $vstrands[$ih]{col};
    ($sum % 2 == 0) ? return 1 : return -1;
}

sub get_direction_stap {
    my $ih = shift;
    my $sum = $vstrands[$ih]{row} + $vstrands[$ih]{col};
    ($sum % 2 == 0) ? return -1 : return 1;
}

## must be same as get_direction_stap()
sub is_upstream {
    my $ih = shift;

    for my $ib (0 .. $#{$vstrands[$ih]->{stap}}) {
	next if $vstrands[$ih]->{skip}[$ib] == -1;
	my $nih = $vstrands[$ih]{stap}[$ib][2];
	my $nib = $vstrands[$ih]{stap}[$ib][3];
	next if ($nib == -1);
	
	if ($ih == $nih) {
	    ($nib > $ib) ? return 1: return -1;
	}
    }
    print STDERR "ERROR: Can't determine direction!\n";
}

sub print_help {
    print STDERR "\nSYNOPSIS\n";
    print STDERR "\tcadnano2pdb.pl --ff=(CHARMM|AMBER) --lattice=(honeycomb|square) [--scaf=SCAFFOLD_SEQUENCE] [--stap=STAPLE_CSV_FILE] [--help] JSON_FILE\n";
    print STDERR "\nPREREQUSITES\n";
    print STDERR "\tThis utility is written in PERL language and requires JSON perl package.\n\n";
    print STDERR "\nDESCRIPTION\n";
    print STDERR "\tThe cadnano2pdb utility converts a JSON file from caDNAno program to a PDB file and write the PDB file to STDOUT. ";
    print STDERR "JSON_FILE is mandatory, and the other options are optional. The options not specified in the commandline will be determined interactively.\n";
    print STDERR "\tScaffold and staple sequence information can be provided using --scaf and --staple options, respectively. ";
    print STDERR "When only one of those two are given, the sequence of the other will be determined complementarily. ";
    print STDERR "\n\tThe options are as follows\n\n";
    print STDERR "\t--help: shows this help information.\n\n";
    print STDERR "\t--ff=(CHARMM|AMBER): determine either CHARMM- and AMBER-compatible PDB outputs. For example, nucleotide and atoms names and the order of atoms in a nucleotide will be consistent with the chosen force field\n\n";
    print STDERR "\t--lattice=(honeycomb|square): specify the lattice type in which JSON file was created.\n\n";
    print STDERR "\t--scaf=SCAFFOLD_SEQUENCE: the sequence of scaffold. The sequence data is taken from the caDNAno program. When not given, staple sequences must be given using --stap option.\n";
    print STDERR "\t\tSCAFFOLD_SEQUENCE=\n";
    for my $s (keys %scafseq) {
	printf STDERR "\t\t\t %-15s (%5d nts)\n", $s, $#{$scafseq{$s}}+1;
    }
    print STDERR "\t--stap=STAPLE_CSV_FILE: CSV file from caDNAno.\n";
}
