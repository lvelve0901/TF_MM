#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;


my $pdb=$ARGV[0]; # the pdb file we are measuring
my $first=$ARGV[1]; # the first strand, input as (chnltr).(resnum):(resnum)  ex: A.1:8
my $second=$ARGV[2]; # the second strand, input as (chnltr).(resnum):(resnum)  ex: A.20:28
my @sstr=($ARGV[3],$ARGV[4]);


my ($rarray,$refmol)=readPDB("iBformDNA.pdb");#iBformDNA_mod.pdb # reference PDB file, needs to be aligned along Z-direction of the frame
my ($carray,$cmol)=readPDB($pdb);

my @chns=(substr($first,0,1),substr($second,0,1));
my @r1=split(/:/,substr($first,2));
my @r2=split(/:/,substr($second,2));

my @res1=($r1[0]);
my @res2=($r2[0]);


# This identifies the residues from the order of the residues in the pdb file, thus it should account 
# for non standard residue numberings that have either letters in the numbers or are non-continuous
my $true=0;
my $iter=0;

my @par1=($sstr[0]=~/([\(]+)([\.]+||\s)([\(]+)/);
my @par2=($sstr[1]=~/([\)]+)([\.]+||\s)([\)]+)/);
if ($par1[1]=~/\s/) {$par1[1]='';}
if ($par2[1]=~/\s/) {$par2[1]='';}
#print length($par1[0])," ",length($par1[1])," ",length($par1[2]),"\n";
#print length($par2[0])," ",length($par2[1])," ",length($par2[2]),"\n";

foreach my $j (@{$carray->{$chns[0]}}) {
  if (!$true) {
    if ($j eq $r1[0]) {$true=1;$iter+=1;
    } else {next;}
  } else {
    $iter+=1;
    if ($iter<=length($par1[0])) {push(@res1,$j);
    } elsif ($iter>(length($par1[0])+length($par1[1])) && $iter<length($par1[0])+length($par1[1])+length($par1[2])) {
      push(@res1,$j);
    } elsif ($iter==(length($par1[0])+length($par1[1])+length($par1[2]))) {push(@res1,$j); last;}
  }
}
$iter=0;
$true=0;
foreach my $j (@{$carray->{$chns[1]}}) {
  if (!$true) {
    if ($j eq $r2[0]) {$true=1;$iter+=1;
    } else {next;}
  } else {
    #print $j,"\n";
    $iter+=1;
    if ($iter<=length($par2[0])) {push(@res2,$j);
    } elsif ($iter>(length($par2[0])+length($par2[1])) && $iter<(length($par2[0])+length($par2[1])+length($par2[2]))) {push(@res2,$j);
    } elsif ($iter==(length($par2[0])+length($par2[1])+length($par2[2]))) {push(@res2,$j); last;}
  }
}
    
if ($res1[-1] ne $r1[1] || $res2[-1] ne $r2[1]) {
  print "@res1    @res2\n";
  die "Residue Count is Inconsistent!!! strand1: $res1[-1],$r1[1] strand2: $res2[-1],$r2[1]\n";}

my $Helix1='';
my $Helix2='';

my $order=1;

if ((length($par1[0])+length($par1[1])+length($par1[2]))<(length($par2[0])+length($par2[1])+length($par2[2]))) {
  # If the first strand is shorter than the second, they are in the wrong order
  $order=0;
} elsif ((length($par1[0])+length($par1[1])+length($par1[2]))==(length($par2[0])+length($par2[1])+length($par2[2]))) {
  # if the lengths of the strands are equal
  if ($res1[0]=~m/^\d+$/ && $res2[0]=~m/^\d+$/) {
    if ($res1[0]>$res2[0]) {
      # if residue numbers are fully numeric and the first residue number of the 
      # first strand is great than the number of the second strand, wrong order
      $order=0;
    } elsif ($res1[0] eq $res2[0] && $chns[0] gt $chns[1]) {
      $order=0;
    }
  } else {
    # residue numbers are not fully numeric
    my ($npt1) = ($res1[0] =~ /(\d+)/);
    my ($npt2) = ($res2[0] =~ /(\d+)/);
    my ($apt1) = ($res1[0] =~ /(\D)/);
    my ($apt2) = ($res2[0] =~ /(\D)/);	
    if ($npt1 > $npt2) {
      $order=0;
    } elsif ($npt1==$npt2) {
      if ($apt1 gt $apt2) {
	$order=0;
      } elsif ($res1[0] eq $res2[0] && $chns[0] gt $chns[1]) {
	$order=0;
      }
    } 
  }
}

if ($order) {
  # they are already in the right order, i.e. first strand is the longer strand  
  # or the length of the two strands is equal and the residue number of the first strand is less than the second
  # or if the residue numbers are the same than the chain letter of the first strand is less than the chain letter of the second
  foreach my $j (0..length($par1[0])-1) {$Helix1=sprintf "%s\_%s.%s,%s.%d",$Helix1,$chns[0],$res1[$j],'Z',10-length($par1[0])+1+$j;}
  foreach my $j ((length($par2[0]))..(length($par2[0])+length($par2[2])-1)) {$Helix1=sprintf "%s\_%s.%s,%s.%d",$Helix1,$chns[1],$res2[$j],'Z',31-length($par2[0])+$j;}  
  foreach my $j ((length($par1[0]))..(length($par1[0])+length($par1[2])-1)) {$Helix2=sprintf "%s\_%s.%s,%s.%d",$Helix2,$chns[0],$res1[$j],'Z',10-length($par1[0])+1+$j;}
  foreach my $j (0..length($par2[0])-1) {$Helix2=sprintf "%s\_%s.%s,%s.%d",$Helix2,$chns[1],$res2[$j],'Z',31-length($par2[0])+$j;}
  #print substr($Helix1,1),"\n";
  #print substr($Helix2,1),"\n";
  #print "@res1\n";
  #print "@res2\n";
  printf "%s %s.%s\:%s,%s.%s\:%s %s.%s\:%s,%s.%s\:%s %s %s %s",$pdb,$chns[0],$res1[0],$res1[length($par1[0])-1],$chns[1],$res2[length($par2[0])],$res2[-1],
    $chns[0],$res1[length($par1[0])],$res1[-1],$chns[1],$res2[0],$res2[length($par2[0])-1],$sstr[0],$sstr[1],abg_measure($cmol,$refmol,substr($Helix1,1),substr($Helix2,1),'back');

} else {
  # we need to reverse them so the longer strand starts on the left
  foreach my $j (0..length($par2[0])-1) {$Helix1=sprintf "%s\_%s.%s,%s.%d",$Helix1,$chns[1],$res2[$j],'Z',10-length($par2[0])+1+$j;}
  foreach my $j ((length($par1[0]))..(length($par1[0])+length($par1[2])-1)) {$Helix1=sprintf "%s\_%s.%s,%s.%d",$Helix1,$chns[0],$res1[$j],'Z',31-length($par1[0])+$j;}  
  foreach my $j ((length($par2[0]))..(length($par2[0])+length($par2[2])-1)) {$Helix2=sprintf "%s\_%s.%s,%s.%d",$Helix2,$chns[1],$res2[$j],'Z',10-length($par2[0])+1+$j;}
  foreach my $j (0..length($par1[0])-1) {$Helix2=sprintf "%s\_%s.%s,%s.%d",$Helix2,$chns[0],$res1[$j],'Z',31-length($par1[0])+$j;} 
  #print $Helix1,"\n";
  printf "%s %s.%s\:%s,%s.%s\:%s %s.%s\:%s,%s.%s\:%s %s %s %s",$pdb,$chns[1],$res2[0],$res2[length($par2[0])-1],$chns[0],$res1[length($par1[0])],$res1[-1],
    $chns[1],$res2[length($par2[0])],$res2[-1],$chns[0],$res1[0],$res1[length($par1[0])-1],$sstr[1],$sstr[0],abg_measure($cmol,$refmol,substr($Helix1,1),substr($Helix2,1),'back');
}


sub readPDB {
  
  my $handle=shift;
  open(PDB,$handle);
  my $residues={};
  my $rlist={};
  my ($lab,$resnumb,$rtype,$atmtype,$chn,$chemtype,$atype);
  my $iter=0;
  my ($anum,$rnum);

 LINE: while (<PDB>) {
    my $ln={};
    ($lab=substr($_,0,6))=~s/ //g;
    last LINE if ($lab=~/END/);

    last LINE if ($lab=~/ENDMDL/);
    next unless ($lab eq "ATOM");
    
   
    ($rtype=substr($_,17,3))=~s/\s+//g;
    if ($rtype eq 'GUA') {$rtype='G';
    } elsif ($rtype eq 'CYT') {$rtype='C';
    } elsif ($rtype eq 'ADE') {$rtype='A';
    } elsif ($rtype eq 'URA') {$rtype='U';}
    
    if (!($rtype eq 'A' || $rtype eq 'G' || $rtype eq 'C' || $rtype eq 'U')) {
	#print STDERR "ERROR: RESTYPE NOT U,A,G,C\n";
	#print STDERR "rtype = $rtype\n";
	#print STDERR "        $_\n";
	next;
    }

    ($chn=substr($_,21,1))=~s/\s+/_/g;  # replace empty chain with _
    ($atmtype=substr($_,12,4))=~s/\s+//g;
    # only take the A model of the atom if there are multiple models...
    if (substr($atmtype,-1) eq 'A') {$atmtype=substr($atmtype,0,-1);($atmtype=substr($_,12,4))=~s/\s+//g;}
    elsif (substr($atmtype,-1) eq 'B') {next;}
    $atmtype=~tr/*/'/;
    ($chemtype=substr($_,72))=~s/[ \n]//g;
    ($anum=substr($_,6,5))=~s/\s+//g;
    ($rnum=substr($_,22,5))=~s/\s+//g;
    
    
    $resnumb=$chn.'.'.$rnum;
    if (!exists($residues->{$resnumb})) {
      $residues->{$resnumb}={};
    }
    
    $ln->{atnum}=$anum;
    $ln->{atype}=$chemtype;
    
    $ln->{xcoor}=substr($_,28,10)+0.0;
    $ln->{ycoor}=substr($_,38,8)+0.0;
    $ln->{zcoor}=substr($_,46,8)+0.0;
    
    $residues->{$resnumb}->{$atmtype}=$ln;
    if (!exists($residues->{$resnumb}->{rtype})) {$residues->{$resnumb}->{rtype}=$rtype;}
    if (!exists($residues->{$resnumb}->{chnltr})) {$residues->{$resnumb}->{chnltr}=$chn;}
    
    if (!exists($rlist->{$chn})) {$rlist->{$chn}=[$rnum];
    } elsif (exists($rlist->{$chn}) && $rlist->{$chn}->[-1] ne $rnum) {
      #if ($chn eq 'A') {print $rnum," ",$rlist->{$chn}->[-1],"\n";}
      push(@{$rlist->{$chn}},$rnum);} 

  }
  close PDB;

  return $rlist,$residues;
}

sub bead_conversion {
  # convert a pdb to the 3-bead representation
  
  my $rna=shift;
  my @reses=keys %$rna;
  my (@sug,@base,@phos,$resnm);
  my ($x,$y,$z);
  my $bead={};
  foreach my $resnum (@reses) {
    
    # the phosphate is the same
    if (!exists($rna->{$resnum}->{P})) {@phos=('NaN','NaN','NaN');
    } else {@phos=($rna->{$resnum}->{P}->{xcoor},$rna->{$resnum}->{P}->{ycoor},$rna->{$resnum}->{P}->{zcoor});}
    
    ($x,$y,$z)=(0,0,0);
    my $num=0;
    # now get the sugar coordinates
    foreach("C1'","C2'","C3'","C4'","O4'","C5'") {
      if (exists($rna->{$resnum}->{$_})) {
	$num+=1;
	$x+=$rna->{$resnum}->{$_}->{xcoor};
	$y+=$rna->{$resnum}->{$_}->{ycoor};
	$z+=$rna->{$resnum}->{$_}->{zcoor};
      }
    }
    if ($num == 6) {@sug=($x/6,$y/6,$z/6);
    } else {@sug=('NaN','NaN','NaN');}
    

    ($x,$y,$z)=(0,0,0);
    $num=0;
    foreach('N1','C2','N3','C4','C5','C6','N7','C8','N9') {
      if (exists($rna->{$resnum}->{$_})) {
	$num+=1;
	$x+=$rna->{$resnum}->{$_}->{xcoor};
	$y+=$rna->{$resnum}->{$_}->{ycoor};
	$z+=$rna->{$resnum}->{$_}->{zcoor};
      }
    }
    if ($num==6 || $num==9) {@base=($x/$num,$y/$num,$z/$num);
    } else {@base=('NaN','NaN','NaN');}
    
    # add to the bead now
    if (!($phos[0] eq 'NaN' && $base[0] eq 'NaN' && $sug[0] eq 'NaN')) {
      $bead->{$resnum}={};
      $bead->{$resnum}->{B}->{xcoor}=$base[0];
      $bead->{$resnum}->{B}->{ycoor}=$base[1];
      $bead->{$resnum}->{B}->{zcoor}=$base[2];

      $bead->{$resnum}->{P}->{xcoor}=$phos[0];
      $bead->{$resnum}->{P}->{ycoor}=$phos[1];
      $bead->{$resnum}->{P}->{zcoor}=$phos[2];
      
      $bead->{$resnum}->{S}->{xcoor}=$sug[0];
      $bead->{$resnum}->{S}->{ycoor}=$sug[1];
      $bead->{$resnum}->{S}->{zcoor}=$sug[2];

      $bead->{$resnum}->{type}=$rna->{$resnum}->{rtype};
      $bead->{$resnum}->{chnltr}=$rna->{$resnum}->{chnltr};
    }
  }
  return $bead;
  
}


sub abg_measure {
    
    # this method performs the same calculation as abg_measureFull,
    # but instead speeds up the calculate by only performing rotations
    # on the sets of atoms of interest

    #my $cmppdb=shift;
    my $compreslist=shift;
    my $refreslist=shift;
    my $h1=shift;
    my $h2=shift;
    my $selemode=shift;

    my $output;
    
    if (!defined($selemode)) {$selemode='back';}

    my $atlist1=&_pairAtoms($compreslist,$refreslist,$h1,$selemode);
    
    if (!defined($atlist1)) {
      die "No defined selection\n";
      return;}

    my @trvector1=&_translateRotate($atlist1);
    if ($trvector1[3]>2) {print STDERR "WARNING! Helix 1 RMSD > 2.0\n";}
    $output=sprintf("HelixRMSD= %.2f",$trvector1[3]);

    my $atlist2=&_pairAtoms($compreslist,$refreslist,$h2,$selemode);
    if (!defined($atlist2)) {
      die "No defined selection\n";
      return;}

    my @trvector2=&_translateRotate($atlist2);
    if ($trvector2[3]>2) {print STDERR "WARNING! Helix 2 RMSD > 2.0\n";}
    $output=sprintf("%s HelixRMSD= %.2f",$output,$trvector2[3]);

    if (!@trvector1 || !@trvector2) {
      die "No defined rotation\n";
      return;}

    my $atlistfinal=[];
    
    for (my $ik=0; $ik<=$#{$atlist2}; $ik++) {

	my $ox=$atlist2->[$ik]->{comp}->{xcoor};
	my $oy=$atlist2->[$ik]->{comp}->{ycoor};
	my $oz=$atlist2->[$ik]->{comp}->{zcoor};

	my $h1x=$ox+$trvector1[0]->[0];
	my $h1y=$oy+$trvector1[0]->[1];
	my $h1z=$oz+$trvector1[0]->[2];

	my $tx1=$trvector1[2]->[1]->[1]*$h1x+$trvector1[2]->[1]->[2]*$h1y+$trvector1[2]->[1]->[3]*$h1z;
	my $ty1=$trvector1[2]->[2]->[1]*$h1x+$trvector1[2]->[2]->[2]*$h1y+$trvector1[2]->[2]->[3]*$h1z;
	my $tz1=$trvector1[2]->[3]->[1]*$h1x+$trvector1[2]->[3]->[2]*$h1y+$trvector1[2]->[3]->[3]*$h1z;
	
	# the atoms that we pass to _translateRotate in the atlist don't need to have
	# all the atom information -- just pass the coordinates...
	my $h1atom={};
	$h1atom->{xcoor}=$tx1+$trvector1[1]->[0];
	$h1atom->{ycoor}=$ty1+$trvector1[1]->[1];
	$h1atom->{zcoor}=$tz1+$trvector1[1]->[2];

	my $h2x=$ox+$trvector2[0]->[0];
	my $h2y=$oy+$trvector2[0]->[1];
	my $h2z=$oz+$trvector2[0]->[2];

	my $tx2=$trvector2[2]->[1]->[1]*$h2x+$trvector2[2]->[1]->[2]*$h2y+$trvector2[2]->[1]->[3]*$h2z;
	my $ty2=$trvector2[2]->[2]->[1]*$h2x+$trvector2[2]->[2]->[2]*$h2y+$trvector2[2]->[2]->[3]*$h2z;
	my $tz2=$trvector2[2]->[3]->[1]*$h2x+$trvector2[2]->[3]->[2]*$h2y+$trvector2[2]->[3]->[3]*$h2z;
	
	my $h2atom={};
	$h2atom->{xcoor}=$tx2+$trvector2[1]->[0];
	$h2atom->{ycoor}=$ty2+$trvector2[1]->[1];
	$h2atom->{zcoor}=$tz2+$trvector2[1]->[2];
	push(@{$atlistfinal},{"ref"=>$h1atom,"comp"=>$h2atom});
    }
    
    my @trvectfinal=&_translateRotate($atlistfinal);
    if ($trvectfinal[3]>=0.01) {
      die "FinalRMSD violation!\n";
      return;
    }

    my $angs=&_rot2eul($trvectfinal[2]);
    #$output=sprintf("%s FinalRMSD= %.2f",$output,$trvectfinal[3]);
    my $minang=&abgminimize($angs->[0],$angs->[1],$angs->[2]);
    $output=sprintf("%s ANGLES:  %.3f  %.3f  %.3f\n",$output,$minang->[0],$minang->[1],$minang->[2]);
    return $output;
}

sub abgminimize {
    # minimize the alpha beta and gamma angles in the Euler space
    my ($a,$b,$g) = @_;
    if($a < -180) { $a+=360; }
    elsif($a > 180) { $a-=360; }
    if($b > 180) { $b-=360; }
    elsif($b < -180) { $b+=360; }
    if($g < -180) { $g+=360; }
    elsif($g > 180) { $g-=360; }
    my ($abg,$d) = ([],[]);
    
    $abg->[0] = [$a,$b,$g];
    $abg->[1] = [$a-180,$b*-1,$g+180];
    $abg->[2] = [$a+180,$b*-1,$g-180];
    $abg->[3] = [$a+180,$b*-1,$g+180];
    $abg->[4] = [$a-180,$b*-1,$g-180];
    
    $d->[0] = sqrt($a**2+$b**2+$g**2);
    $d->[1] = sqrt(($a-180)**2+($b*-1)**2+($g+180)**2);
    $d->[2] = sqrt(($a+180)**2+($b*-1)**2+($g-180)**2);
    $d->[3] = sqrt(($a+180)**2+($b*-1)**2+($g+180)**2);
    $d->[4] = sqrt(($a-180)**2+($b*-1)**2+($g-180)**2);
    
    my ($pck,$i,$j) = (1e+10,0,0);
    foreach (@$d) { if($pck > $_) { $pck = $_; $i = $j; }; $j++; }
    return ($abg->[$i]);
}

sub _rot2eul {
    # Deconvolute the rotation matrix to the alpha, beta, gamma angles
    my $rot = shift;
    my $eul = [];
    my $pi=4*atan2(1,1);
    $eul->[0] = (180.0/$pi)*atan2($rot->[3]->[2],$rot->[3]->[1]);
    $eul->[2] = (180.0/$pi)*atan2($rot->[2]->[3],$rot->[1]->[3]*-1.0);
    my $tmp = $rot->[2]->[3]/sin($eul->[2]*($pi/180.0));
    $eul->[1] = (180.0/$pi)*atan2($tmp,$rot->[3]->[3]);

    return $eul;
}

sub _pairAtoms {

    ## $paired is an array that can take on various formats, though it each 
    ## case it requires the specification of two residues that we wish to match
    ## Selection mode specifices the types of atoms within each residue
    ## that we want to match.
    ## $paired format is taken as a string argument.  It can either be an
    ## individual list of paired residues -- pairs separated by commas, sets 
    ## separated by _
    ##  (strnum).(resnum),(strnum).(resnum)_(strnum).(resnum),(strnum).(resnum)...
    ## $paired can also be input as string of continious residues, with multiple
    ## pairings also strung together by _:
    ##  (strnum).(resnum_start):(resnum_end),(strnum).(resnum_start):(resnum_end)_...
    ## When listed this way the numberings must be continous, and the matches will
    ## be made such that (resnum_start) is matched to (resnum_start), and moving along
    ## the residues
    ## NOTE: By convention the reference residues are listed secondly

    my $cmpres=shift;
    my $refres=shift;
    my $paired=shift;
    my $selmode=shift;
    
    my $atmlist=();
    
    $selmode="all" if (!defined $selmode || $selmode eq "");
    $selmode=lc $selmode;

    my $search;
    $search="^(P|S)\$" if ($selmode eq "threebead" || $selmode eq 'toy');
    $search="^P\$"      if ($selmode eq "p");
    $search="^(B|S)\$"      if ($selmode eq "bs");
    $search="^.*\$"      if ($selmode eq "all" || $selmode eq "heavy");
    $search="^(P|O3'|C3'|C4'|C5'|O5'|C2'|C1'|O4'|O2')\$" if ($selmode eq "back");
    

    # parse the $paired argument
    my $rlist=[];
    my %firsts;
    my @selects=split(/_/,$paired);
    my $counter=0;
    foreach my $i (@selects) {
      #print $i,"\n";
	my @tmp=split(/,/,$i);

	if ($tmp[0]=~/:/) {
	    $tmp[0] =~/^([A-Za-z0-9]).(\d+):(\d+)/;
	    my $str1=$1;
	    my @bnd1=($2,$3);
	    $firsts{$str1.'.'.$bnd1[0]}='';
	    $tmp[1] =~/^([A-Za-z0-9]).(\d+):(\d+)/;
	    my $str2=$1;
	    my @bnd2=($2,$3);
	    $firsts{$str2.'.'.$bnd2[0]}='';

	    if (abs($bnd2[0]-$bnd2[1]) != abs($bnd1[0]-$bnd1[1])) {
	      #print "Selections are not the same length!!!\n";
	      die "Selections are not the same length!!!\n"; 
	      next;}
	    my $b1inc=(($bnd1[1]-$bnd1[0])>0)?1:-1;
	    my $b2inc=(($bnd2[1]-$bnd2[0])>0)?1:-1;
	    foreach (0..abs($bnd1[0]-$bnd1[1])) {
		push(@{$rlist},[$str1.'.'.($bnd1[0]+$b1inc*$_),$str2.'.'.($bnd2[0]+$b2inc*$_)]);
	    }
	} else {
	  push(@{$rlist},[$tmp[0],$tmp[1]]);
	  if ($counter==0 || $counter==3) {
	    $firsts{$tmp[0]}='';
	    $firsts{$tmp[1]}='';
	  }
	}
	$counter+=1;
    }

    # set up the $atmlist hash
    for (my $i=0; $i<=$#{$rlist};$i++) {
	if (!exists($cmpres->{$rlist->[$i]->[0]}) || !exists($refres->{$rlist->[$i]->[1]})) {
	    print $rlist->[$i]->[0].".  .".$rlist->[$i]->[1],"\n";
	    die "Selected Residue does not exist!!!". "\n";
	}
	   
	foreach (keys %{$cmpres->{$rlist->[$i]->[0]}}) {
	    if ($_=~/$search/) {
		if (!exists($refres->{$rlist->[$i]->[1]}->{$_})) {
		    printf "cannot match atom $_ in residue ".$rlist->[$i]->[1]."\n";
		} elsif ($_ eq 'P' && exists($firsts{$rlist->[$i]->[0]}) && $selmode eq 'toy') {next;
		} else {
		    push(@{$atmlist},{'comp'=>$cmpres->{$rlist->[$i]->[0]}->{$_},
				      'ref'=>$refres->{$rlist->[$i]->[1]}->{$_}});
		}
	    }
	}
    }
    return $atmlist;
}

sub _translateRotate {
    
    ## This method is modified from the _translateRoate method of the
    ## Analyze.pm package of the MMTSB toolset.  
    ## Specifically, this method is modified to allow fitting of atoms
    ## /residues that do not have the same numbering
    
    ## Note: $keylist is pre-paired by like atoms, so selection is done
    ## in _pairAtoms

    my $pairlist=shift;
    if ($#{$pairlist}<0) {die "no atoms selected\n";}

    my ($cxref,$cyref,$czref);
    my ($cxcmp,$cycmp,$czcmp);
    $cxref=$cyref=$czref=$cxcmp=$cycmp=$czcmp=0.0;

    # $pairlist comes in as an hashref of hashes with the paired residue pointers
    # $pairlist=[
    #           {'comp'=>cmpatm0,'ref'=>refatm0},
    #           {'comp'=>cmpatm1,'ref'=>refatm1},
    #           etc... }

    my $nref=0;
    for (my $ik=0; $ik<=$#{$pairlist}; $ik++) {
	my $rat=$pairlist->[$ik]->{ref};
	my $cat=$pairlist->[$ik]->{comp};
	$cxref+=$rat->{xcoor};
	$cyref+=$rat->{ycoor};
	$czref+=$rat->{zcoor};
	
	$cxcmp+=$cat->{xcoor};
	$cycmp+=$cat->{ycoor};
	$czcmp+=$cat->{zcoor};
	
	$nref++;
    }

    $cxref/=$nref;
    $cyref/=$nref;
    $czref/=$nref;
    
    $cxcmp/=$nref;
    $cycmp/=$nref;
    $czcmp/=$nref;

    # init the matrix
    my $r=();
    for (my $i=1; $i<=3; $i++) {
	$r->[$i]=();
	for (my $j=1; $j<=3; $j++) {
	    $r->[$i]->[$j]=0.0;
	}
    }
    
    for (my $ik=0; $ik<=$#{$pairlist}; $ik++) {
	my $cat=$pairlist->[$ik]->{comp};
	my $rat=$pairlist->[$ik]->{ref};
	my $cx=$cat->{xcoor}-$cxcmp;
	my $cy=$cat->{ycoor}-$cycmp;
	my $cz=$cat->{zcoor}-$czcmp;
	my $rx=$rat->{xcoor}-$cxref;
	my $ry=$rat->{ycoor}-$cyref;
	my $rz=$rat->{zcoor}-$czref;
	
	
	$r->[1]->[1]+=$cx*$rx;
	$r->[2]->[1]+=$cx*$ry;
	$r->[3]->[1]+=$cx*$rz;
	
	$r->[1]->[2]+=$cy*$rx;
	$r->[2]->[2]+=$cy*$ry;
	$r->[3]->[2]+=$cy*$rz;
	
	$r->[1]->[3]+=$cz*$rx;
	$r->[2]->[3]+=$cz*$ry;
	$r->[3]->[3]+=$cz*$rz;
    }
    
    my $u=_frotu($r);
    my $translate1=[-1*$cxcmp,-1*$cycmp,-1*$czcmp];
    my $translate2=[$cxref,$cyref,$czref];



    # rotate the compare structure!
    my $rms=0;
    my $numdiff=0;
    for (my $ik=0; $ik<=$#{$pairlist}; $ik++) {
	my $cat=$pairlist->[$ik]->{comp};
	my $rat=$pairlist->[$ik]->{ref};
	
	my $x=$cat->{xcoor}+$translate1->[0];
	my $y=$cat->{ycoor}+$translate1->[1];
	my $z=$cat->{zcoor}+$translate1->[2];

	my $tx=$u->[1]->[1]*$x+$u->[1]->[2]*$y+$u->[1]->[3]*$z;
	my $ty=$u->[2]->[1]*$x+$u->[2]->[2]*$y+$u->[2]->[3]*$z;
	my $tz=$u->[3]->[1]*$x+$u->[3]->[2]*$y+$u->[3]->[3]*$z;

	my $dx=$rat->{xcoor}-($tx+$translate2->[0]);
	my $dy=$rat->{ycoor}-($ty+$translate2->[1]);
	my $dz=$rat->{zcoor}-($tz+$translate2->[2]);
	$rms+=$dx*$dx+$dy*$dy+$dz*$dz;
	$numdiff+=1;
    }
    
    
    return ($translate1,$translate2,$u,sqrt($rms/$numdiff));
  
}

sub _frotu {
  ## This method is taken wholesale from the MMTSB Analyze.pm package

  my $r=shift;

  my ($i,$j,$k);

  my $det=0.0;
  for ($i=1; $i<=3; $i++) {
    my $i1=$i+1;
    $i1-=3 if ($i1>3);
    my $i2=$i+2;
    $i2-=3 if ($i2>3);
    $det+=$r->[$i]->[1]*($r->[$i1]->[2]*$r->[$i2]->[3]-$r->[$i2]->[2]*$r->[$i1]->[3]);
  }

  my $ipt=0;
  my $w;

  for ($i=1; $i<=3; $i++) {
    for ($j=$i; $j<=3; $j++) {
      $ipt++;
      $w->[$ipt]=0.0;
      for ($k=1; $k<=3; $k++) {
	$w->[$ipt]+=$r->[$j]->[$k]*$r->[$i]->[$k];
      }
    }
  }

  my $trace=$w->[1]+$w->[4]+$w->[6];

  my $u=();

  if ($trace<3.0E-6) {
    for ($i=1; $i<=3; $i++) {
      $u->[$i]=();
      for ($j=1; $j<=3; $j++) {
	$u->[$i]->[$j]=0.0;
      }
      $u->[$i]->[$i]=1.0;
    }
    return $u;
  }

  my ($vec,$ev)=_diagq(3,3,$w);

  my $a=();
  $ipt=1;
  for ($i=1; $i<=3; $i++) {
    for ($j=1; $j<=3; $j++) {
      $a->[$j]->[$i]=$vec->[$ipt];
      $ipt++;
    }
  }

  for ($i=1; $i<=3; $i++) {
    $ev->[$i]=sqrt(abs($ev->[$i]));
    $ev->[$i]=1.0E-6 if ($ev->[$i]<1.0E-6);
  }

  $ev->[1]=-$ev->[1] if ($det<0.0);

  my $b=();
  for ($j=1; $j<=3; $j++) {
    $b->[$j]=();
  }

  for ($j=1; $j<=3; $j++) {
    my $evs=$ev->[$j];
    for ($i=1; $i<=3; $i++) {
      $b->[$i]->[$j]=0.0;
      for ($k=1; $k<=3; $k++) {
	$b->[$i]->[$j]+=$r->[$k]->[$i]*$a->[$k]->[$j]/$evs;
      }
    }
  }

  $det=0.0;
  for ($i=1; $i<=3; $i++) {
    my $i1=$i+1;
    $i1-=3 if ($i1>3);
    my $i2=$i+2;
    $i2-=3 if ($i2>3);
    $det+=$a->[$i]->[1]*($a->[$i1]->[2]*$a->[$i2]->[3]-$a->[$i2]->[2]*$a->[$i1]->[3]);
  }

  for ($j=1; $j<=3; $j++) {
    if (abs($ev->[$j]) <= 1.0E-6) {
      my $jp=$j+1;
      my $jq=$j+2;
      $jp-=3 if ($jp>3);
      $jq-=3 if ($jq>3);
      for ($k=1; $k<=3; $k++) {
	my $kp=$k+1;
	my $kq=$k+2;
	$kp-=3 if ($kp>3);
	$kq-=3 if ($kq>3);
	$b->[$k]->[$j]=$b->[$kp]->[$jp]*$b->[$kq]->[$jq]-$b->[$kp]->[$jq]*$b->[$kq]->[$jp];
	$b->[$k]->[$j]=-$b->[$k]->[$j] if ($det<0.0);
      }
    }

    my $c=0.0;
    for ($k=1; $k<=3; $k++) {
      $c+=$b->[$k]->[$j]*$b->[$k]->[$j];
    }

    $c=($c>1.0E-10)?1.0/sqrt($c):0.0;

    for ($k=1; $k<=3; $k++) {
      $b->[$k]->[$j]*=$c;
    }
  }

  for ($j=1; $j<=3; $j++) {
    $u->[$i]=();
    for ($i=1; $i<=3; $i++) {
      $u->[$i]->[$j]=0.0;
      for ($k=1; $k<=3; $k++) {
	$u->[$i]->[$j]+=$a->[$i]->[$k]*$b->[$j]->[$k];
      }
    }
  }

  for ($j=1; $j<=3; $j++) {
    my $c=0.0;
    for ($k=1; $k<=3; $k++) {
      $c+=$u->[$k]->[$j]*$u->[$k]->[$j];
    }
    $c=($c>1.0E-10)?1.0/sqrt($c):0.0;

    for ($k=1; $k<=3; $k++) {
      $u->[$k]->[$j]*=$c;
    }
  }

  $det=0.0;
  
  for ($i=1; $i<=3; $i++) {
    my $i1=$i+1;
    $i1-=3 if ($i1>3);
    my $i2=$i+2;
    $i2-=3 if ($i2>3);
    $det+=$u->[$i]->[1]*($u->[$i1]->[2]*$u->[$i2]->[3]-$u->[$i2]->[2]*$u->[$i1]->[3]);
  }
  
  if (abs($det-1.0)>1.0E-4) {
    printf "non-unitary rotation matrix, determinant: %f\n",$det;
    return;
  }
  else {return $u;}

}    

sub _diagq {
  my $nx=shift;
  my $nfrqx=shift;
  my $dd=shift;

  my $vec=();
  my $ev=();

  my @a;
  my @b;
  my @p;
  my @w;
  my @ta;
  my @tb;
  my @y;
  
  my $nadd=0;

  my $eta=2.22045E-16;
  my $theta=4.4923E+307;

  my $n=$nx;
  my $nev=$nfrqx;
  my $nevadd=$nev+$nadd; 

  my $del1=$eta/100.0;
  my $delta=$eta*$eta*100.0;
  my $small=$eta*$eta/100.0;
  my $delbig=$theta*$delta/1000.0;
  my $theta1=1000.0/$theta;
  my $toler=100.0*$eta;
  my $rpower=8388608.0;
  my $rpow1=$rpower*0.50;
  my $rand1=$rpower-3.0;
  my $dunity=1.0;

  my $factor=0.0;
  my $ntot=int(($n*($n+1))/2);

  my ($i,$j,$k,$l,$m);

  for ($i=1; $i<=$ntot; $i++) {
    $factor=abs($dd->[$i]) if ($factor<abs($dd->[$i]));
  }

  if ($factor<$theta1) {
    printf "zero matrix passed in _diagq\n";

    for ($i=1; $i<=$nev; $i++) {
      $ev->[$i]=0.0;
      my $ipt=($i-1)*$n;
      for ($j=1; $j<$n; $j++) {
	$ipt++;
	$vec->[$ipt]=($i+$nadd == $j)?1.0:0.0;
      }
    }
    return ($vec,$ev);
  }

  my $ij=0;
  my $anorm=0.0;
  
  for ($i=1; $i<=$n; $i++) {
    for ($j=$i; $j<=$n; $j++) {
      $ij++;
      my $u=($dd->[$ij]/$factor)*($dd->[$ij]/$factor);
      $u*=0.5 if ($i == $j);
      $anorm+=$u;
    }
  }
  $anorm=sqrt($anorm+$anorm)*$factor;
  my $anormr=$dunity/$anorm;

  for ($i=1; $i<=$ntot; $i++) {
    $dd->[$i]*=$anormr;
  }

  my $nn=$n-1;
  my $mi=0;
  my $mi1=$n-1;

  for ($i=1; $i<=$nn; $i++) {
    my $sum1=0.0;
    $b[$i]=0.0;
    my $ji=$i+1;
    my $ipt=$mi+$i;
    $a[$i]=$dd->[$ipt];
    $ipt++;
    my $bx=$dd->[$ipt];
    my $ji2=$ji+1;
    if ($ji == $n) {
      $b[$i]=$bx;
      $dd->[$mi+$ji]=0.0;
      $mi+=$mi1;
      $mi1--;
    } else {
      for ($j=$ji2; $j<=$n; $j++) {
	$ipt++;
	$sum1+=$dd->[$ipt]*$dd->[$ipt];
      }
      
      if ($sum1<=$small) {
	$b[$i]=$bx;
	$dd->[$mi+$ji]=0.0;
	$mi+=$mi1;
	$mi1--;
      } else {
	my $s=sqrt($sum1+$bx*$bx);
	my $sgn=($bx>=0.0)?abs($dunity):-abs($dunity);
	my $temp=abs($bx);
	$w[$ji]=sqrt(0.5*($dunity+($temp/$s)));
	$ipt=$mi+$ji;
	$dd->[$ipt]=$w[$ji];
	my $ii=$i+2;
	if ($ii<=$n) {
	  $temp=$sgn/(2.0*$w[$ji]*$s);
	  for ($j=$ii; $j<=$n; $j++) {
	    $ipt++;
	    $w[$j]=$temp*$dd->[$ipt];
	    $dd->[$ipt]=$w[$j];
	  }
	}

	$b[$i]=-$sgn*$s;
	
	for ($j=$ji; $j<=$n; $j++) {
	  $p[$j]=0.0;
	}

	my $ml=$mi+$mi1;
	my $ml1=$mi1-1;


	for ($l=$ji; $l<=$n; $l++) {
	  $ipt=$ml+$l;
	  for ($m=$l; $m<=$n; $m++) {
	    $bx=$dd->[$ipt];
	    $p[$l]+=$bx*$w[$m];
	    $p[$m]+=$bx*$w[$l] if ($l!=$m);
	    $ipt++;
	  }
	  $ml+=$ml1;
	  $ml1--;
	}
	my $xkap=0.0;
	
	for ($k=$ji; $k<=$n; $k++) {
	  $xkap+=$w[$k]*$p[$k];
	}

	for ($l=$ji; $l<=$n; $l++) {
	  $p[$l]-=$xkap*$w[$l];
	}
	
	my $mj=$mi+$mi1;
	my $mj1=$mi1-1;

	for ($j=$ji; $j<=$n; $j++) {
	  for ($k=$j; $k<=$n; $k++) {
	    my $expr=($p[$j]*$w[$k])+($p[$k]*$w[$j]);
	    $dd->[$mj+$k]-=$expr+$expr;
	  }
	  $mj+=$mj1;
	  $mj1--;
	}
	
	$mi+=$mi1;
	$mi1--;
      }
    }
  }

  $a[$n]=$dd->[$mi+$n];
  $b[$n]=0.0;

  my $alimit=1.0;
  for ($i=1; $i<=$n; $i++) {
    $w[$i]=$b[$i];
    $b[$i]*=$b[$i];
  }
  
  for ($i=1; $i<=$nevadd; $i++) {
    $ev->[$i]=$alimit;
  }
  my $rootl=-$alimit;

  for ($i=1; $i<=$nevadd; $i++) {  
    my $rootx=$alimit;
    for ($j=$i; $j<=$nevadd; $j++) {
      $rootx=$ev->[$j] if ($ev->[$j]<$rootx);
    }
    $ev->[$i]=$rootx;

    my $trial=($rootl+$ev->[$i])*0.5;

    while(abs($trial-$rootl) > 1.0E-15 && abs($trial-$ev->[$i]) > 1.0E-15) {
      my $nomtch=$n;
      $j=1;
      
      do {
	my $f0=$a[$j]-$trial;
	#print 

	while ($j<=$n && abs($f0)>=$theta1) {
	  $nomtch-- if ($f0>=0.0);
	  $j++;
          if($j>$n) { next } 
	  # note: this was an error from the MMTSB toolset where $j=4 and $a[4] is undefined
	  # throwing a warning. However it should be noted that this doesn't matter.  
	  # The $f0 quantity is only used in the context of executing this while statement
	  # Thus, on the final increment of j, where j=4 and then $j>$n so the while loop will
	  # exit anyways, so evaluation of $f0 is not needed. 
	  $f0=$a[$j]-$trial-$b[$j-1]/$f0;
	}  
	if ($j<=$n) {
	  $j+=2;
	  $nomtch--;
	}
      }	while ($j<=$n);

      if ($nomtch>=$i) {
	$ev->[$i]=$trial;
	my $nom=($nevadd<=$nomtch)?$nevadd:$nomtch;
	$ev->[$nom]=$trial;
      } else {
	$rootl=$trial;
      }

      $trial=($rootl+$ev->[$i])*0.5;      
    }
  }

  for ($i=1; $i<=$nev; $i++) {
    $ev->[$i]=$ev->[$i+$nadd];
  }

  my $ia=0;
  
  for ($i=1; $i<=$nev; $i++) {
    my $aroot=$ev->[$i];
    for ($j=1; $j<=$n; $j++) {
      $y[$j]=1.0;
    }
    $ia=-1 if ($i==1 || abs($ev->[$i-1]-$aroot)>=$toler);
    $ia++;

    my $elim1=$a[1]-$aroot;
    my $elim2=$w[1];

    for ($j=1; $j<=$nn; $j++) {
      my $temp;
      if (abs($elim1)<=abs($w[$j])) {
	$ta[$j]=$w[$j];
	$tb[$j]=$a[$j+1]-$aroot;
	$p[$j]=$w[$j+1];
	$temp=(abs($w[$j])>$theta1)?$elim1/$w[$j]:1.0;
	$elim1=$elim2-$temp*$tb[$j];
	$elim2=-$temp*$w[$j+1];
      } else {
	$ta[$j]=$elim1;
	$tb[$j]=$elim2;
	$p[$j]=0.0;
	$temp=$w[$j]/$elim1;
	$elim1=$a[$j+1]-$aroot-$temp*$elim2;
	$elim2=$w[$j+1];
      }
      $b[$j]=$temp;
    }

    $ta[$n]=$elim1;
    $tb[$n]=0.0;
    $p[$n]=0.0;
    $p[$nn]=0.0;
    my $iter=1;
    
    if ($ia!=0) { 
      for ($j=1; $j<=$n; $j++) {
	my $rand1=(4099.0*$rand1 % $rpower);
	$y[$j]=$rand1/$rpow1-1.0;
      }
    }

    do {
      $l=$n+1;
      
      for ($j=1; $j<=$n; $j++) {
	$l--;
	do {
	  if (($n-$l-1)<0) {
	    $elim1=$y[$l];
	  } elsif (($n-$l-1)==0) {
	    $elim1=$y[$l]-$y[$l+1]*$tb[$l];
	  } else {
	    $elim1=$y[$l]-$y[$l+1]*$tb[$l]-$y[$l+2]*$p[$l];
	  }
	  if ($elim1>$delbig || $elim1<-$delbig) {
	    for ($k=1; $k<=$n; $k++) {
	      $y[$k]=$y[$k]/$delbig;
	    }
	  } 	  
	} while ($elim1>$delbig || $elim1<-$delbig);
	
	my $temp=$ta[$l];
	$temp=$delta if (abs($temp)<$delta);
	$y[$l]=$elim1/$temp;
      }
      
      if ($iter==1) {
	$elim1=$y[1];
	for ($j=1; $j<=$nn; $j++) {
	  if (abs($ta[$j]-$w[$j])<1E-15) {
	    $y[$j]=$y[$j+1];
	    $elim1=$elim1-$y[$j+1]*$b[$j];
	  } else {
	    $y[$j]=$elim1;
	    $elim1=$y[$j+1]-$elim1*$b[$j];
	  }
	}
	$y[$n]=$elim1;
      } 
      $iter++;
    } while ($iter<=2);
    
    my $ipt;
    if ($ia != 0) {
      for (my $j1=1; $j1<=$ia; $j1++) {
	$k=$i-$j1;
	my $temp=0.0;
	$ipt=($k-1)*$n;
	for ($j=1; $j<=$n; $j++) {
	  $ipt++;
	  $temp+=$y[$j]*$vec->[$ipt];
	}
	$ipt=($k-1)*$n;
	for ($j=1; $j<=$n; $j++) {
	  $ipt++;
	  $y[$j]-=$temp*$vec->[$ipt];
	}
      }
    }
    
    $elim1=0.0;
    for ($j=1; $j<=$n; $j++) {
      $elim1=abs($y[$j]) if (abs($y[$j])>$elim1);
    }
    my $temp=0.0;
    for ($j=1; $j<=$n; $j++) {
      $elim2=$y[$j]/$elim1;
      $temp+=$elim2*$elim2;
    }
    $temp=$dunity/(sqrt($temp)*$elim1);
    for ($j=1; $j<=$n; $j++) {
      $y[$j]*=$temp;
      $y[$j]=0.0 if (abs($y[$j])<$del1);
    }
    $ipt=($i-1)*$n;
    for ($j=1; $j<=$n; $j++) {
      $ipt++;
      $vec->[$ipt]=$y[$j];
    }
  }

  my $ipt;
  my $kk;
  for ($i=1; $i<=$nev; $i++) {
    $ipt=($i-1)*$n;
    for ($j=1; $j<=$n; $j++) {
      $ipt++;
      $y[$j]=$vec->[$ipt];
    }
    $l=$n-2;
    my $mk=($n*($n-1))/2-3;
    my $mk1=3;
    
    my $t;
    for ($j=1; $j<=$l; $j++) {
      $t=0.0;
      $k=$n-$j-1;
      $m=$k+1;
      for ($kk=$m; $kk<=$n; $kk++) {
	$t+=$dd->[$mk+$kk]*$y[$kk];
      }
      for ($kk=$m; $kk<=$n; $kk++) {
	my $epr=$t*$dd->[$mk+$kk];
	$y[$kk]-=$epr+$epr;
      }
      $mk-=$mk1;
      $mk1++;
    }

    $t=0.0;
    for ($j=1; $j<=$n; $j++) {
      $t+=$y[$j]*$y[$j];
    }
    
    my $xnorm=sqrt($t);
    my $xnorm1=$dunity/$xnorm;
    for ($j=1; $j<=$n; $j++) {
      $y[$j]*=$xnorm1;
    }
    
    $ipt=($i-1)*$n;
    for ($j=1; $j<=$n; $j++) {
      $ipt++;
      $vec->[$ipt]=$y[$j];
    }
  }

  for ($i=1; $i<=$n; $i++) {
    $ev->[$i]=$ev->[$i]*$anorm;
  }

  return ($vec,$ev);
}
