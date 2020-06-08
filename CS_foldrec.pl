#!/usr/bin/perl
#===============================================================
# FOLD RECOGNITION :: COMPLEMENTARITY SCORE
# CS_gl, CS_cp 
# Computed from side-chain shape complementarity (Sm) & electrostatic complementarity (Em) of buried amino acids at the globular protein interior 
# Burial Cutoff: 0.30
# Sequence -> Threaded onto a fold -> (Sm, Em) -> CSgl, CScp
# ---------------------------------------------------------------
# Reference: 
# ARTICLE| VOLUME 102, ISSUE 11, P2605-2614, JUNE 06, 2012
# Self-Complementarity within Proteins: Bridging the Gap between Binding and Folding
# Sankar Basu, Dhananjay Bhattacharyya, Rahul Banerjee
# Open ArchiveDOI:https://doi.org/10.1016/j.bpj.2012.04.029
# --------------------------------------------------------------
#===============================================================
#


$ref="
# ---------------------------------------------------------------
# Reference: 
# ARTICLE| VOLUME 102, ISSUE 11, P2605-2614, JUNE 06, 2012
# Self-Complementarity within Proteins: Bridging the Gap between Binding and Folding
# Sankar Basu, Dhananjay Bhattacharyya, Rahul Banerjee
# Open ArchiveDOI:https://doi.org/10.1016/j.bpj.2012.04.029
# --------------------------------------------------------------
#===============================================================\n";


# .Sm
$ifile1 = $ARGV[0] || die "--------------------------------------------------------\nEnter Sm profile (\$basename.Sm) & Em profile (\$basename.Em) as arg-1 & arg-2;\nUsage: ./CS_foldrec.pl <\$basename>.Sm <\$basename>.Em\n--------------------------------------------------------\n";	
chomp $ifile1;

# .Em
$ifile2 = $ARGV[1] || die "--------------------------------------------------------\nEnter Sm profile (\$basename.Sm) & Em profile (\$basename.Em) as arg-1 & arg-2;\nUsage: ./CS_foldrec.pl <\$basename>.Sm <\$basename>.Em\n--------------------------------------------------------\n";
chomp $ifile2;

$burcutoff = 0.30;

$tot1 = 141954;
$tot2 = 84288;
$totsq = $tot1*$tot2;

#========================== libraries ===================

open (LIB1,"<Scomp.mlib");
open (LIB2,"<Ecomp.mlib");

open (RGB1,"<resGbur1.prob");
open (RGB2,"<resGbur2.prob");
open (RGB3,"<resGbur3.prob");
open (RGB4,"<resGbur4.prob");

open (BGR1,"<burGres1.prob");
open (BGR2,"<burGres2.prob");
open (BGR3,"<burGres3.prob");
open (BGR4,"<burGres4.prob");

open (SBR1,"<flcb.sctab");
open (SBR2,"<flpb.sctab");
open (SBR3,"<flpe.sctab");
open (SBR4,"<flfe.sctab");

open (EBR1,"<flcb.ectab");
open (EBR2,"<flpb.ectab");
open (EBR3,"<flpe.ectab");
open (EBR4,"<flfe.ectab");

open (ABR,"<burAVG.mlib");

open (LTZ1,"<lorentz1.mlib");
open (LTZ2,"<lorentz2.mlib");
open (LTZ3,"<lorentz3.mlib");
open (LTZ4,"<lorentz4.mlib");

@ltz1 = <LTZ1>;
@ltz2 = <LTZ2>;
@ltz3 = <LTZ3>;
@ltz4 = <LTZ4>;

@rgb1 = <RGB1>;
@rgb2 = <RGB2>;
@rgb3 = <RGB3>;
@rgb4 = <RGB4>;

@bgr1 = <BGR1>;
@bgr2 = <BGR2>;
@bgr3 = <BGR3>;
@bgr4 = <BGR4>;

@sclib = <LIB1>;
@eclib = <LIB2>;

@sbr1 = <SBR1>;
@sbr2 = <SBR2>;
@sbr3 = <SBR3>;
@sbr4 = <SBR4>;

@ebr1 = <EBR1>;
@ebr2 = <EBR2>;
@ebr3 = <EBR3>;
@ebr4 = <EBR4>;

@abr = <ABR>;

foreach (@sclib) {chomp $_;}
foreach (@eclib) {chomp $_;}

foreach (@rgb1) {chomp $_;}
foreach (@rgb2) {chomp $_;}
foreach (@rgb3) {chomp $_;}
foreach (@rgb4) {chomp $_;}

foreach (@bgr1) {chomp $_;}
foreach (@bgr2) {chomp $_;}
foreach (@bgr3) {chomp $_;}
foreach (@bgr4) {chomp $_;}

foreach (@sbr1) {chomp $_;}
foreach (@sbr2) {chomp $_;}
foreach (@sbr3) {chomp $_;}
foreach (@sbr4) {chomp $_;}

foreach (@ebr1) {chomp $_;}
foreach (@ebr2) {chomp $_;}
foreach (@ebr3) {chomp $_;}
foreach (@ebr4) {chomp $_;}

foreach (@abr) {chomp $_;}

foreach (@ltz1) {chomp $_;}
foreach (@ltz2) {chomp $_;}
foreach (@ltz3) {chomp $_;}

$lsc = @sclib;
$lec = @eclib;

#print "$lsc  $lec\n";

if ($lsc >= $lec)
{
$lh = $lsc;
}
else
{
$lh = $lec;
}

#========================================================

$pi = 4*atan2(1.00,1.00);

$K1 = sqrt(2*$pi);
$K2 = sqrt(2*$pi);

open (INP1,"<$ifile1");
@dat1 = <INP1>;

open (INP2,"<$ifile2");
@dat2 = <INP2>;

$l = @dat1;
$l2 = @dat2;

if ($l != $l2)
{
die "Length missmatch: $ifile1: $l || $ifile2: $l2\n";
}

$cnt = 0;

$Esum = 0.00;
$PsumP = 0.00;
$sum_scaled = 0.00;
$sumSclSig = 0.00;

$Essum = 0.0;
$Eesum = 0.0;
$Lssum = 0.0;
$Lesum = 0.0;
$Lsum = 0.0;
$GLsum = 0.0;

$Psumbur = 0.00;
$Psumbgr = 0.00;

$Psumbsc = 0.00;
$Psumbec = 0.00;

$sumN = 0.00;
$sumTOT = 0.00;

$pnS1 = 0;
$pnS2 = 0;
$pnS3 = 0;

$pnE1 = 0;
$pnE2 = 0;
$pnE3 = 0;

for $i (0..$l-1)
{
chomp $dat1[$i];
chomp $dat2[$i];

$rn1 = int(substr($dat1[$i],0,3));
$res1 = substr($dat1[$i],4,3);
$bur1 = substr($dat1[$i],9,4);
$scsc = substr($dat1[$i],25,6);

$rn2 = int(substr($dat2[$i],0,3));
$res2 = substr($dat2[$i],4,3);
$bur2 = substr($dat2[$i],12,4);
$scec = substr($dat2[$i],31,6);

$btag = '';
$pfilesc = '';
$pfileec = '';

#print "$scsc  $scec\n";

$Esc = 0.00;		# Empiriacl Gaussian type
$Eec = 0.00;

$ps = 0.00;		# P(Sc|Br,Res)
$pe = 0.00;		# P(Ec|Br,Res)

$pbur = 0.00;		# P(Res|Br)

$pbursc = 0.00;		# P(Res|{Br,Sc})
$pburec = 0.00;		# P(Res|{Br,Ec})

$nts = 0.00;
$nte = 0.00;
$ns = 0.00;
$ne = 0.00;
$nt = 0.00;
$num = 0.00;
$p = 0.00;


	if ($bur1 <= 0.05)
	{
	$btag = 'cb';
	$pfilesc = $res1.'1'.'.sctab';
	$pfileec = $res2.'1'.'.ectab';
	}
	elsif ($bur1 > 0.05 && $bur1 <= 0.15)
	{
	$btag = 'pb';
	$pfilesc = $res1.'2'.'.sctab';
	$pfileec = $res2.'2'.'.ectab';
	}
	elsif ($bur1 > 0.15 && $bur1 <= 0.30)
	{
	$btag = 'pe';
	$pfilesc = $res1.'3'.'.sctab';
	$pfileec = $res2.'3'.'.ectab';
	}
	elsif ($bur1 > 0.30)
	{
	$btag = 'fe';
	}

	if ($bur2 <= 0.05)
	{
	$btag2 = 'cb';
	}
	elsif ($bur2 > 0.05 && $bur2 <= 0.15)
	{
	$btag2 = 'pb';
	}
	elsif ($bur2 > 0.15 && $bur2 <= 0.30)
	{
	$btag2 = 'pe';
	}
	elsif ($bur2 > 0.30)
	{
	$btag2 = 'fe';
	}

#print $btag,"\n";

	if ($res1 eq $res2 && $rn1 == $rn2 && $res1 ne 'GLY')
	{
		if ($btag ne $btag2)
		{
		die "Missmatch in burial status: $ifile1: $rn1  $res1  $bur1 || $ifile2: $rn2  $res2  $bur2\n";		
		}
		if ($bur1 <= $burcutoff && $bur2 <= $burcutoff)
		{
		$cnt++;
			foreach $bn (@abr)
			{
			chomp $bn;
				if (substr($bn,0,3) eq $res1)
				{
				$mbur = substr($bn,7,7);
				}
			}			
		$ScBur = abs($bur1-$mbur)/$mbur;		# |Br(i)-<Br>|/<Br>

		open (PSC,"<$pfilesc");		# related Sc bin file 
		open (PEC,"<$pfileec");		# related Ec bin file
		@psc = <PSC>;
		@pec = <PEC>;
		close PSC;
		close PEC;
#===========================================================================================
#=================================== P(Res|Br), P(Res|{Br,Sc}), P(Res|{Br,Ec})  ============
#===========================================================================================
			if ($btag eq 'cb')
			{
			# Lorentz parameter
				foreach $lz (@ltz1)
				{
					if ($res1 eq substr($lz,0,3))
					{
					$Emed = substr($lz,8,5);		# Median
					$HW = substr($lz,18,5);			# Half width at half maximum hight
#					print "$res1  $btag  $Emed  $HW\n";
					}
				}
				foreach $x (@rgb1)
				{
					if (substr($x,0,3) eq $res1)
					{
					$pbur = substr($x,16,6);		# P(Res|Br)
					}
				}
				foreach $x (@bgr1)
				{
					if (substr($x,0,3) eq $res1)
					{
					$pbr = substr($x,22,5);	# P(Br|Res)
					}
				}
				foreach $y (@sbr1)
				{
				$lbt = substr($y,0,8);
				$ubt = substr($y,10,8);				
				$nt = substr($y,19,6);					
					if ($scsc > $lbt && $scsc <= $ubt)
					{
					$nts = $nt;
					}
				}
				foreach $z (@ebr1)
				{
				$lbt = substr($z,0,8);
				$ubt = substr($z,10,8);				
				$nt = substr($z,19,6);
					if ($scec > $lbt && $scec <= $ubt)
					{
					$nte = $nt;
					}
				}
			}
			elsif ($btag eq 'pb')
			{
			# Lorentz parameter
				foreach $lz (@ltz2)
				{
					if ($res1 eq substr($lz,0,3))
					{
					$Emed = substr($lz,8,5);		# median 
					$HW = substr($lz,18,5);			# Half width at half maximum hight
#					print "$res1  $btag  $Emed  $HW\n";
					}
				}
				foreach $x (@rgb2)
				{
					if (substr($x,0,3) eq $res1)
					{
					$pbur = substr($x,16,6);		# P(Res|Br)
					}
				}
				foreach $x (@bgr2)
				{
					if (substr($x,0,3) eq $res1)
					{
					$pbr = substr($x,22,5);	# P(Br|Res)
					}
				}
				foreach $y (@sbr2)
				{
				$lbt = substr($y,0,8);
				$ubt = substr($y,10,8);				
				$nt = substr($y,19,6);
					if ($scsc > $lbt && $scsc <= $ubt)
					{
					$nts = $nt;
					}
				}
				foreach $z (@ebr2)
				{
				$lbt = substr($z,0,8);
				$ubt = substr($z,10,8);				
				$nt = substr($z,19,6);
					if ($scec > $lbt && $scec <= $ubt)
					{
					$nte = $nt;
					}
				}
			}
			elsif ($btag eq 'pe')
			{
			# Lorentz parameter
				foreach $lz (@ltz3)
				{
					if ($res1 eq substr($lz,0,3))
					{
					$Emed = substr($lz,8,5);		# Median
					$HW = substr($lz,18,5);			# Half width at half maximum hight
#					print "$res1  $btag  $Emed  $HW\n";
					}
				}
				foreach $x (@rgb3)
				{
					if (substr($x,0,3) eq $res1)
					{
					$pbur = substr($x,16,6);		# P(Res|Br)
					}
				}
				foreach $x (@bgr3)
				{
					if (substr($x,0,3) eq $res1)
					{
					$pbr = substr($x,22,5);	# P(Br|Res)
					}
				}
				foreach $y (@sbr3)
				{
				$lbt = substr($y,0,8);
				$ubt = substr($y,10,8);				
				$nt = substr($y,19,6);
					if ($scsc > $lbt && $scsc <= $ubt)
					{
					$nts = $nt;
					}
				}
				foreach $z (@ebr3)
				{
				$lbt = substr($z,0,8);
				$ubt = substr($z,10,8);				
				$nt = substr($z,19,6);
					if ($scec > $lbt && $scec <= $ubt)
					{
					$nte = $nt;
					}
				}
			}
			elsif ($btag eq 'fe')
			{
			# Lorentz parameter
				foreach $lz (@ltz4)
				{
					if ($res1 eq substr($lz,0,3))
					{
					$Emed = substr($lz,8,5);		# median 
					$HW = substr($lz,18,5);			# Half width at half maximum hight
#					print "$res1  $btag  $Emed  $HW\n";
					}
				}
				foreach $x (@rgb4)
				{
					if (substr($x,0,3) eq $res1)
					{
					$pbur = substr($x,16,6);		# P(Res|Br)
					}
				}
				foreach $x (@bgr4)
				{
					if (substr($x,0,3) eq $res1)
					{
					$pbr = substr($x,22,5);	# P(Br|Res)
					}
				}
				foreach $y (@sbr4)
				{
				$lbt = substr($y,0,8);
				$ubt = substr($y,10,8);				
				$nt = substr($y,19,6);
					if ($scsc > $lbt && $scsc <= $ubt)
					{
					$nts = $nt;
					}
				}
				foreach $z (@ebr4)
				{
				$lbt = substr($z,0,8);
				$ubt = substr($z,10,8);				
				$nt = substr($z,19,6);
					if ($scec > $lbt && $scec <= $ubt)
					{
					$nte = $nt;
					}
				}
			}
#================================================================================
#============================== P(Sc|{Br,Res}) ==================================
#============================== P(Ec|{Br,Res}) ==================================
#================================================================================
			foreach $s (@psc)
			{
			chomp $s;
			$lb = substr($s,0,8);
			$ub = substr($s,10,8);
			$p = substr($s,30,5);
			$num = substr($s,19,6);
				if ($scsc == -1.00)	# minimum possible value for sc
				{
				$ps = substr($psc[0],30,5);
				}
				else
				{
					if ($scsc > $lb && $scsc <= $ub)
					{
					$ps = $p;
					$ns = $num;
					}
				}
			}

			foreach $e (@pec)
			{
			chomp $e;
			$lb = substr($e,0,8);
			$ub = substr($e,10,8);
			$p = substr($e,30,5);
			$num = substr($e,19,6);
				if ($scec == -1.00)	# minimum possible value for ec
				{
				$pe = substr($pec[0],30,5);
				}
				else
				{
					if ($scec > $lb && $scec <= $ub)				
					{
					$pe = $p;
					$ne = $num;
					}
				}
			}
#			print "$rn1  $res1     $btag     $scsc     $scec     $ns  $nts  $ps   $ne  $nte  $pe  $Emed  $HW\n";
#================================================================================
#========================= Read from main Library ===============================
#================================================================================
			for $n (0..$lh-1)
			{
			$rt1 = substr($sclib[$n],0,3);
			$br1 = substr($sclib[$n],8,5);
			$Scm = substr($sclib[$n],39,6);
			$sig1 = substr($sclib[$n],48,5);
			$bind1 = '';

			$rt2 = substr($eclib[$n],0,3);
			$br2 = substr($eclib[$n],8,5);
			$Ecm = substr($eclib[$n],39,6);
			$sig2 = substr($eclib[$n],48,5);
			$bind2 = '';

				if ($br1 <= 0.05)
				{
				$bind1 = 'cb';
				}
				elsif ($br1 > 0.05 && $br1 <= 0.15)
				{
				$bind1 = 'pb';
				}
				elsif ($br1 > 0.15 && $br1 <= 0.30)
				{
				$bind1 = 'pe';			
				}
				elsif ($br1 > 0.30)
				{
				$bind1 = 'fe';
				}
#			print $bind1,"\n";

				if ($br2 <= 0.05)
				{
				$bind2 = 'cb';
				}
				elsif ($br2 > 0.05 && $br2 <= 0.15)
				{
				$bind2 = 'pb';
				}
				elsif ($br2 > 0.15 && $br2 <= 0.30)
				{
				$bind2 = 'pe';
				}
				elsif ($br2 > 0.30)
				{
				$bind2 = 'fe';
				}
#			print $bind2,"\n";

#================================================================================
#=================== empirical deviation maximizing measures ====================
#================================================================================

				if ($res1 eq $rt1 && $btag eq $bind1)
				{
				$Esc = exp(-0.5*(($scsc-$Scm)/$sig1)**2)/($K1*$sig1);		# Gaussian
				$Scaled_Sc = abs(($scsc-$Scm)/$Scm);
				$ScaleSigSc = abs(($scsc-$Scm)/$sig1);
#============================================================================================
#  Penalty Scheme for SC
#============================================================================================
				$lcS = $Scm - 2*$sig1;
				$ucS = $Scm + 2*$sig1;				
					if ($scsc < $lcS || $scsc > $ucS)
					{
					#print "$lcS  $ucS  $scsc\n";
						if ($btag eq 'cb')
						{
						$pnS1++;
						}
						elsif ($btag eq 'pb')
						{
						$pnS2++;
						}
						elsif ($btag eq 'pe')
						{
						$pnS3++;
						}
					}						
				}

				if ($res1 eq $rt2 && $btag eq $bind2)
				{
				$Eec = exp(-0.5*(($scec-$Ecm)/$sig2)**2)/($K2*$sig2);		# Gaussian
				$Scaled_Ec = abs(($scec-$Ecm)/$Ecm);
				$ScaleSigEc = abs(($scec-$Ecm)/$sig2);
				$Lec = ($HW/(($HW)**2 + ($scec-$Emed)**2))/$pi;			# Lorentzian
#============================================================================================
#  Penalty Scheme for EC
#============================================================================================
				$lcE = $Ecm - $sig2;
				$ucE = $Ecm + $sig2;				
					if ($scec < 0)		# For anticomplementarity
					{
#					print "$lcE  $ucE  $scec\n";
						if ($btag eq 'cb')
						{
						$pnE1++;
						}
						elsif ($btag eq 'pb')
						{
						$pnE2++;
						}
						elsif ($btag eq 'pe')
						{
						$pnE3++;
						}
					}
				}
			}
		}
	}

	if ($nts > 0)
	{
	$pbursc = $ns/$nts;		# P(Res|Br,Sc)
	}
	else 
	{
	$pbursc = 0.00;
	}

	if ($nte > 0)
	{
	$pburec = $ne/$nte;		# P(Res|Br,Ec)
	}
	else
	{
	$pburec = 0.00;
	}

#================== Gaussian ====================
	
$Essum = $Essum + $Esc;
$Eesum = $Eesum + $Esc;	
	
$Eprod = $Esc*$Eec;
$Esum = $Esum + $Eprod;

#================== Lorentzian ==================

$Lesum = $Lesum + $Lec;

#================== Gaussian (Sc) * Lorentzian (Ec) ==

$GLprod = $Esc*$Lec;
$GLsum = $GLsum + $GLprod;

#======================Cond Prob======================

$psumSc = $psumSc + $ps;
$psumEc = $psumEc + $pe;

$prodP = $ps*$pe;
$PsumP = $PsumP + $prodP;

#============================================

$Psumbur = $Psumbur + $pbur;
$Psumbgr = $Psumbgr + $pbr;

$Psumbsc = $Psumbsc + $pbursc;

$prodPR = $pbursc*$pburec;

$Psumbscec = $Psumbscec + $prodPR;

$sumN = $sumN + ($ns*$ne);

$prod_scaled = $Scaled_Sc*$Scaled_Ec;			# * $ScBur;
$sum_scaled = $sum_scaled + $prod_scaled;

$pSclSig = $ScaleSigSc*$ScaleSigEc;
$sumSclSig = $sumSclSig + $pSclSig;

$r1 = $ns/$tot1;		# total probability (P(res|br,Sc)) = N(res|br,Sc)/N(tot1)
$r2 = $ne/$tot2;		# total probability (P(res|br,Ec)) = N(res|br,Ec)/N(tot2)

$prodr1r2 = $r1*$r2;
$sumTOT = $sumTOT + $prodr1r2;

#printf "%3d  %3s  %2s  %6.3f  %6.3f  %5d  %5d  %8.5f  %5d  %5d  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",$rn1,$res1,$btag,$scsc,$scec,$ns,$nts,$pbursc,$ne,$nte,$pburec,$prodPR,$ps,$pe,$prodP,$Esc,$Eec,$Eprod;
#printf "%3d  %6.3f  %6.3f  %5d  %5d  %8.5f  %5d  %5d  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",$rn1,$scsc,$scec,$ns,$nts,$pbursc,$ne,$nte,$pburec,$prodPR,$ps,$pe,$prodP,$Esc,$Eec,$Eprod;

	if ($bur1 <= $burcutoff && $res1 ne 'GLY')
	{
#	printf "%3d %3s %2s %6.3f %6.3f %5d %5d %6d %5d %5d %8.5f %8.5f %8.5f %8.5f %8.5f %9.7f %6.3f %6.3f\n",$rn1,$res1,$btag,$scsc,$scec,$ns,$ne,sqrt($totsq),$nts,$nte,$pbursc,$pburec,$prodPR,$r1,$r2,$prodr1r2,$ps,$pe;
	}
}

#print "$Psumbur  $cnt\n";

$En = $Esum/$cnt;
$PnP = sprintf("%8.5f",($PsumP/$cnt));
$PnP=~s/\s//g;

$PnB = $Psumbur/$cnt;
$PnBgR = $Psumbgr/$cnt;

$Pnbs = $Psumbsc/$cnt;
$Pnbse = $Psumbscec/$cnt;

$raw = $sumN/$cnt;
$nraw = $raw/$totsq;

$NormTOT = $sumTOT/$cnt;

$pS1 = $psumSc/$cnt;
$pE2 = $psumEc/$cnt;

$Gs = $Essum/$cnt;
$Ge = $Eesum/$cnt;

$Le = $Lesum/$cnt;

$GsLe = sprintf("%8.3f",($GLsum/$cnt));
$GsLe=~s/\s//g;
#========== ADD PENALTY =======================

#printf "%8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %3d\n",$Gs,$Ge,$En,$Le,$GsLe,$pS1,$pE2,$PnP,$cnt;
print "-----------------------------------\n";
print "-----------------------------------\n\n";
print "CSgl=$GsLe  CScp=$PnP  Nbur=$cnt\n\n";
print "-----------------------------------\n";
print "-----------------------------------\n";

print "-----------------------------------------------\n";
print "Average values reported in proteins:\n";
print "CSgl: 3.7 (±0.437); CScp: 0.015 (±0.0017)\n";
print "Mean-3*sd: 2.4 (CSgl); 0.01 (CScp) (baseline)\n";
print "-----------------------------------------------\n";

print $ref,"\n";



