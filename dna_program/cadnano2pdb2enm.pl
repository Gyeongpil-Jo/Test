#!/usr/bin/perl -w

use Getopt::Long;
use constant PI => 4 * atan2 1, 1;

use strict;
use warnings;
use List::Util 'shuffle';
use Math::Trig;
use Math::Trig ':radial';


my $help;
my $namd;
my $debug;
my $cut = 6.0;
my $k=1;    # kcal/mol/A^2

my @word = (
);
   
sub usage {
    print STDERR @word;
    exit;
}

GetOptions(
    "help" => \$help,
    "namd" => \$namd,
    "k=s" => \$k,
    "cut=s" => \$cut,   # cut in angstrom
    "debug" => \$debug
);

&usage if defined $help;

my %hb_atoms = (
    ADE  => [ 'N6', 'N1' ],  #, "H2"
    THY  => [ 'O4', 'N3' ],  #, "O2"
    GUA  => [ 'O6', 'N1' ],  #, "N2"
    CYT  => [ 'N4', 'N3' ],  #, "O2"
    CYTH => [ 'N4', 'N3' ]   #, "O2"
);

my @pdb;
my $old_ch = "xxx";
my $old_resi = 9999;
my $old_resn = "xxx";
my %data;
my $nn=0;
while (<>) {
    next unless (/^(ATOM  )(.{5}).(.{4})(.)(.{4})(.)(.{4}).{4}(.{8})(.{8})(.{8})(.*)/);
    ++$nn;

    my ($ii,$atom,$resn,$ch,$resi,$x,$y,$z) = ($2,$3,$5,$6,$7,$8,$9,$10);
    #next if ($_ =~ / H/);   # shit. chain id can be H.
    next if ($atom =~ /(P|O1P|O2P|H|')/);
    #next unless ($atom =~ /(N1|C2|N3|C4|C5|C6|N7|C8|N9)/);
    next unless ($atom =~ /(N1|C2|N3|C4|C5|C6|N7|C8|N9|O4|O2|N6|H2|N4|O6|N2)/);
    #print "atom=$atom index=$nn resi=$resi\n";

    $atom =~ s/ //g;
    $resn =~ s/[53 ]//g;

    ## resname
    next unless ($resn =~ /(ADE|DA|THY|DT|GUA|DG|CYT|DC)/);
    if    ($resn =~ /^(ADE|DA|A)/) { $resn = "ADE"; } 
    elsif ($resn =~ /^(THY|DT|T)/) { $resn = "THY"; }
    elsif ($resn =~ /^(GUA|DG|G)/) { $resn = "GUA"; }
    elsif ($resn =~ /^(CYT|DC|C)/) { $resn = "CYT"; }

    if ($resi != $old_resi) {
	#print "new chain $old_resi => $resi\n" if ($resi < $old_resi);
	push @pdb, {%data} if %data;

	%data = ( resn => $resn );
	$old_resi = $resi;
    }

    $data{$atom}{xyz} = [$x,$y,$z];
    $data{$atom}{index} = $nn;
}

## push last data.
push @pdb, {%data} if %data;

printf STDERR "Number of nucleotides = %d\n", $#pdb+1;

### calculate COM of base
for my $ri (0 .. $#pdb) {
    my @com = ();
    my $nnn = 0;
    #for my $ai (keys %{$pdb[$ri]}) {
    for my $ai ( 'N1', 'C2', 'N3', 'C4', 'C5', 'C6' ) {
	next if $ai eq "resn";
	#print "$ch..$ri..$ai\n";
	$com[0] += $pdb[$ri]{$ai}{xyz}[0];
	$com[1] += $pdb[$ri]{$ai}{xyz}[1];
	$com[2] += $pdb[$ri]{$ai}{xyz}[2];
	++$nnn;
    }
    $pdb[$ri]{com} = [$com[0]/$nnn,$com[1]/$nnn,$com[2]/$nnn];
    #printf STDERR "%f,%f,%f\n", $com[0]/$nnn,$com[1]/$nnn,$com[2]/$nnn;
}

my $FENM;
my $FHBN;
if (defined $namd) {
    open $FENM, ">enm.extrabonds";
    open $FHBN, ">hbonds.extrabonds";
}
else {
    open $FENM, ">enm.stack.itp";
    open $FHBN, ">enm.hbonds.itp";
}
### make network
print $FENM "[ bonds ]\n" unless defined $namd;
print $FHBN "[ bonds ]\n" unless defined $namd;
my ($resni, $resnj);
for my $ri (0 .. $#pdb) {
    $resni = $pdb[$ri]{resn};
    my @hbi = @{$hb_atoms{$resni}};
    for my $rj (($ri+1) .. $#pdb) {
	$resnj = $pdb[$rj]{resn};
	my @hbj = @{$hb_atoms{$resnj}};
	my $rd = vec_norm(vec_sub(@{$pdb[$ri]{com}},@{$pdb[$rj]{com}}));
	next if ($rd > 20);
	#printf STDERR "rd = %7.3f $resni($ri) $resnj($rj)\n", $rd;
	#printf STDERR "$resni($pdb[$ri]{'C2'}{index}) = (%f,%f,%f) ", @{$pdb[$ri]{'C2'}{xyz}};
	#printf STDERR "com = (%f,%f,%f)\n", @{$pdb[$ri]{com}};
	#printf STDERR "$resnj($pdb[$rj]{'C2'}{index}) = (%f,%f,%f) ", @{$pdb[$rj]{'C2'}{xyz}};
	#printf STDERR "com = (%f,%f,%f)\n", @{$pdb[$rj]{com}};

	#######################################
	# determine if H-bond pair.
	#print STDERR "$resni, $resnj\n";
	#print STDERR @hbj;
	my $dd_major  = vec_norm(vec_sub(@{$pdb[$ri]{$hbi[0]}{xyz}},
		                         @{$pdb[$rj]{$hbj[0]}{xyz}})      );
	my @hb_vec_center = vec_sub(@{$pdb[$ri]{$hbi[1]}{xyz}},
		                    @{$pdb[$rj]{$hbj[1]}{xyz}});
	my $dd_center = vec_norm(@hb_vec_center);

	my $fh = $FENM;
	my $myk = $k;	# k for ENM.
	my $mycut = $cut;	# cut for ENM.

	if ( ($dd_major < 3.5) && ($dd_center < 3.5) ) {
	    # potential H-bond pair. need to examine further.
	    #print STDERR "$pdb[$ri]{'C2'}{index}, $pdb[$ri]{'C5'}{index}\n";
	    #print STDERR "$pdb[$rj]{'C2'}{index}, $pdb[$rj]{'C5'}{index}\n";

	    # see if base planes are in plane.
	    # compute C5-C2-C2-C5 dihedral angle.
	    #my @c2c5 = vec_sub(@{$pdb[$ri]{"C2"}{xyz}}, @{$pdb[$ri]{"C5"}{xyz}});
	    #my @c5c5 = vec_sub(@{$pdb[$rj]{"C2"}{xyz}}, @{$pdb[$ri]{"C5"}{xyz}});
	    #my $phi = do_pdih( @{$pdb[$ri]{'C5'}{xyz}}, @{$pdb[$ri]{'C2'}{xyz}},
	    #	               @{$pdb[$rj]{'C2'}{xyz}}, @{$pdb[$rj]{'C5'}{xyz}} );
	    #printf STDERR "%4d %4d, d= %5.2f phi= %5.2f\n",$ri, $rj, $dd_center, $phi;
	    my @cm2cm = vec_unit(vec_sub(@{$pdb[$ri]{com}}, @{$pdb[$rj]{com}}));
	    # get normal vector of ri
	    my @c2n1 = vec_unit(vec_sub(@{$pdb[$ri]{'C2'}{xyz}}, @{$pdb[$ri]{'N1'}{xyz}}));
	    my @n3c2 = vec_unit(vec_sub(@{$pdb[$ri]{'N3'}{xyz}}, @{$pdb[$ri]{'C2'}{xyz}}));
	    my @nvec = vec_cross(@c2n1, @n3c2);
	    my $inner = vec_inner(@cm2cm, @nvec);
	    #printf STDERR "%4d %4d, d= %5.2f inner= %5.2f\n",$ri, $rj, $dd_center, $inner;
	    #printf STDERR "%4d %4d\n", $ri, $rj if abs($inner) < 0.05;
	    
	    if (abs($inner) < 0.05) {
		$fh = $FHBN;
		$myk = 1;
		$mycut = 12;
		# pair detection is not complete. no difference for now.
		$myk = $k;
		$mycut = $cut;
		print STDERR "Residue $ri and $rj are pair\n";
	    }
	}

	#######################################
	# print ENM or HBN
	for my $ai (keys %{$pdb[$ri]}) {
	    next if $ai eq "resn";
	    next if $ai eq "com";
	    my @xyzi = @{$pdb[$ri]{$ai}{xyz}};
	    for my $aj (keys %{$pdb[$rj]}) {
		next if $aj eq "resn";
		next if $aj eq "com";
		my $dd = vec_norm(vec_sub(@xyzi,@{$pdb[$rj]{$aj}{xyz}}));

		next if $dd > $mycut;

		if (!defined $namd) {	# gromacs
		    printf $fh ("%10d%10d%10d%10.3g%10.3g\n", 
			$pdb[$ri]{$ai}{index}, $pdb[$rj]{$aj}{index},6, $dd*0.1, $myk*418.4);
		}
		else {	# namd: index is zero-based
		    printf $fh ("bond%10d%10d%10.3g%10.3g\n", 
			$pdb[$ri]{$ai}{index}-1, $pdb[$rj]{$aj}{index}-1,$myk,$dd);
		}
	    }
	}
    }
}

close $FHBN;
close $FENM;

exit;

#print "[ moleculetype ]\n";
#printf "%s\t\t%d\n","DNA",1;
#print "\n\n";
#
#print "[ atoms ]\n";
#print ";   nr       type  resnr residue  atom   cgnr     charge       mass\n";
#for my $i (0 .. $#atom) {
#    ### remove spaces
#    $atom[$i] =~ s/ //g;
#    $resn[$i] =~ s/ //g;
#
#    ### atomtype
#    my $type = $rtp{$resn[$i]}{$atom[$i]};
#    unless (defined $type) {
#	print STDERR "atom ..$atom[$i]..$resn[$i].. not found\n";
#    }
#
#    my $chg = ($atom[$i] =~ /^P/) ? -1 : 0;
#
#    # search for the element mass from the name and then
#    # store the string suitable for outputing one line for
#    # this atom type
#    my $mass;
#
#    foreach my $mass_hash ( \%masses, \%singlemasses ) {
#	$mass = find_mass $mass_hash, $atom[$i];
#	last if defined $mass;
#    }
#    die "Didn't find mass for atom type $atom[$i]\n" unless defined $mass;
#
#    printf "%5d %8s %5d %8s %8s %8d %10.3f %10.3f\n", 
#	$i+1,$type,$resi[$i],$resn[$i],$atom[$i],$i+1,$chg,$mass;
#}
#print "\n\n";
#
#print "; ENM with cut = $cut nm\n";
#print "[ bonds ]\n";
#for my $i (0 .. $#atom-1) {
#    for my $j ($i+1 .. $#atom) {
#	my $d = (defined $pbc) ? 
#	    pbc_dx($x[$i],$y[$i],$z[$i],$x[$j],$y[$j],$z[$j]) :
#	    vec_norm(vec_sub($x[$i],$y[$i],$z[$i],$x[$j],$y[$j],$z[$j]));
#
#	$d /= 10;   # nm
#	next if ($d > $cut);
#	printf "%10d%10d%10d%10.3g%10d\n", $i+1, $j+1, 1, $d, 5000;
#    }
#}


exit;



#sub deg2rad { PI * $_[0] / 180.0 }

sub vec_inner {
    return ($_[0]*$_[3]+$_[1]*$_[4]+$_[2]*$_[5]);
}

sub vec_cross {
    return ($_[1]*$_[5]-$_[2]*$_[4],
            $_[2]*$_[3]-$_[0]*$_[5],
            $_[0]*$_[4]-$_[1]*$_[3]);
}

sub vec_sub {
    return ($_[0]-$_[3],$_[1]-$_[4],$_[2]-$_[5]);
}

sub vec_norm {
    return sqrt($_[0]*$_[0]+$_[1]*$_[1]+$_[2]*$_[2]);
}

sub vec_unit {
    my $len = vec_norm(@_);

    return ($_[0]/$len,$_[1]/$len,$_[2]/$len);
}

# gromacs dih_angle()
sub do_pdih {
    my @xi = splice(@_,0,3);
    my @xj = splice(@_,0,3);
    my @xk = splice(@_,0,3);
    my @xl = splice(@_,0,3);

    ## HJ-specific: 
    ## translate j and k so that extend ij & kl vectors by factor of 1.5
    #my @ji = vec_sub(@xj,@xi);
    #@xj = ($xi[0]+$ji[0]*1.5, $xi[1]+$ji[1]*1.5, $xi[2]+$ji[2]*1.5);
    #my @kl = vec_sub(@xk,@xl);
    #@xk = ($xl[0]+$kl[0]*1.5, $xl[1]+$kl[1]*1.5, $xl[2]+$kl[2]*1.5);

    my @r_ij = vec_sub(@xi,@xj);
    my @r_kj = vec_sub(@xk,@xj);
    my @r_kl = vec_sub(@xk,@xl);
    #print STDERR vec_norm(@ji)," ",vec_norm(@r_ij),"\n";

    my @m = vec_cross(@r_ij,@r_kj);
    my @n = vec_cross(@r_kj,@r_kl);
    my @w = vec_cross(@m,@n);
    my $wlen = vec_norm(@w);
    my $s = vec_inner(@m,@n);

    my $phi = rad2deg(atan2($wlen,$s)); ## -180 to +180

    my $ipr = vec_inner(@r_ij,@n);
    my $sign = ($ipr<0.0)?-1:1;

    return $phi*$sign;
}
