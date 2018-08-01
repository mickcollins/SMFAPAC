#!/usr/bin/perl

$llim[1]=$ARGV[0];
$ulim[1]=$ARGV[1];

$n=$ARGV[2];

$CODEDIR=`awk 'NR==1,NR==1 {print \$1}' CODENAME`;
chomp($CODEDIR);
$CODEDIR=~ s/^\s+|\s+$//g;

$EXEDIR=`awk 'NR==1,NR==1 {print \$1}' EXENAME`;
chomp($EXEDIR);
$EXEDIR=~ s/^\s+|\s+$//g;

$deriv =`awk 'NR==2,NR==2 {print \$1}' IN_JOBTYPE`;
chomp($deriv);
$deriv=~ s/^\s+|\s+$//g;

# some data needed in getcoords_qch
@atomic_labels=('H','He','Li','Be','B','C','N','O','F','Ne',
'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb',
'Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',
'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os',
'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
'Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',
'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds');
$numlabels=$#atomic_labels;

if ($llim[1] > 0) {
$junk=`rm -f FragDerivatives$n`;

 for (my $ic=$llim[1]; $ic<=$ulim[1]; $ic++) {
  open(FR,">>FragDerivatives$n");
  my $mf=`wc -l COORD$ic`;
  my @fields = split / /, $mf;
  $natomf=$fields[0];
  $natomf = $natomf - 1; 
  print FR "$natomf\n";
  close FR;
  &getcoords_qch("FRAG$ic.log");
  system("cat thesecoords$n >> FragDerivatives$n");
  &getenergy_qch("FRAG$ic.log");
  open(FR,">>FragDerivatives$n");
  print FR "$TotEn\n";
  close FR;
  if ($deriv >= 1) {
   &getforces_qch("FRAG$ic.log");
   system("cat theseforces$n >> FragDerivatives$n");
  }
  if ($deriv >= 2) {
   &gethessian_qch("FRAG$ic.com.fchk");
   system("cat hessian$n >> FragDerivatives$n");
   &getdipdr_qch("FRAG$ic.log");
   system("echo '$natomf' >> FragDipderivs$n");
   system("cat thesedipdr$n >> FragDipderivs$n");
  }
 }
}


# now the ab initio nonbonded jobs
$AllTot = 0;
if ($deriv >= 1) {
 system("rm -f abforces$n");
}
if ($deriv >= 2) {
 system("rm -f abhessians$n");
}
open($fh,"<abfiles.$n");
while (<$fh>) {
$row=$_;
chomp($row);
 my @fields=split(/\./,$row);
 $ic1=$fields[1];
 $ic2=$fields[2];
 &getenergy_qch("ab.$ic1.$ic2.log");
 $TotE=$TotEn;
# $absign=`awk '/Isg_coeff/ {print \$2}' ab.$ic1.$ic2.com`;
 system("awk /coeff/ ab.$ic1.$ic2.com > frub$n");
 $absign=`awk 'NR==1,NR==1 {print \$2}' frub$n`;
 chomp $absign;
 system("rm -f frub$n");
 if ($deriv >= 1) {
  open(TMP1,">>abforces$n");
  print TMP1 "$ic1","  ","$ic2\n";
  print TMP1 "$absign\n";
  &getforces_qch("ab.$ic1.$ic2.log");
  $nat1=$numberofatoms;
  print TMP1 "$nat1\n";
  close TMP1;
  system("cat theseforces$n >> abforces$n");
 }
 if ($deriv >= 2) {
  open(TMP2,">>abhessians$n");
  print TMP2 "$ic1","  ","$ic2\n";
  print TMP2 "$absign\n";
  &gethessian_qch("ab.$ic1.$ic2.com.fchk",$nat1);
  print TMP2 "$nat1\n";
  close TMP2;
  system("cat hessian$n >> abhessians$n");
 }

RETRY1:
system("mkdir nbdr$ic1 2> /dev/null");
if( $? > 8 ) {
system("sleep 0.028");
goto RETRY1;
}else{
 &getenergy_qch("nb.$ic1.0.log");
 $TotE=$TotE - $TotEn;
 if ($deriv >= 1) {
  &getforces_qch("nb.$ic1.0.log");
  $nat2=$numberofatoms;
  open(TMP1,">>abforces$n");
  print TMP1 "$ic1\n";
  print TMP1 "$nat2\n";
  close TMP1;
  system("cat theseforces$n >> abforces$n");
 }
 if ($deriv >= 2) {
  &gethessian_qch("nb.$ic1.0.com.fchk",$nat2);
  open(TMP2,">>abhessians$n");
  print TMP2 "$ic1\n";
  print TMP2 "$nat2\n";
  close TMP2;
  system("cat hessian$n >> abhessians$n");
 }
 $junk=`rmdir nbdr$ic1 2> /dev/null`;
}

RETRY2:
system("mkdir nbdr$ic2 2> /dev/null");
if( $? > 8 ) {
system("sleep 0.028");
goto RETRY2;
}else{
&getenergy_qch("nb.$ic2.0.log");
 $TotE=$TotE - $TotEn;
 if ($deriv >= 1) {
  &getforces_qch("nb.$ic2.0.log");
  $nat3=$numberofatoms;
  open(TMP1,">>abforces$n");
  print TMP1 "$ic2\n";
  print TMP1 "$nat3\n";
  close TMP1;
  system("cat theseforces$n >> abforces$n");
 }
 if ($deriv >= 2) {
  &gethessian_qch("nb.$ic2.0.com.fchk",$nat3);
  open(TMP2,">>abhessians$n");
  print TMP2 "$ic2\n";
  print TMP2 "$nat3\n";
  close TMP2;
  system("cat hessian$n >> abhessians$n");
 }
$junk=`rmdir nbdr$ic2 2> /dev/null`;
}

$AllTot = $AllTot + $absign * $TotE;

}
close TMP1;
close TMP2;
open (TMP1,">NearInt_energy$n");
print TMP1 "$AllTot\n";
close TMP1;

sub getcoords_qch {
 open(FRAGJOB,"<$_[0]");
 open(TMP,">thesecoords$n");
 my $skip=0;
 my $writit=0;
    $firsttime=0;
  while (<FRAGJOB>) {
   if (/Standard Nuclear Orientation/) {
    $skip=$. + 3;
   }  else {
    if (/Z-Matrix orientation/) {
     $skip=$. + 5;
    }
   }
 if ($. == $skip) {$writit=1};
 if($firsttime == 1) {$writit=0};
 if (/-----------------/) {
   if($writit == 1) {$firsttime=1};
   $writit=0;
 }

 if( $writit == 1) {
  $line1="$_";
  @line=split(/\s+/,$line1);
# get the atomic number for each element
  $numele=0;
  for($ia=0;$ia <= $numlabels;$ia++) {
  if($line[2] eq $atomic_labels[$ia] ) {$numele=$ia+1};
  }
  print TMP "0    $numele    0","   $line[3]","   $line[4]","   $line[5]\n";
      }
  }
 close TMP;
 close FRAGJOB;
}

sub getdipdr_qch {
open(FRAGJOB,"<$_[0]");
open(TMP,">thesedipdr$n");
$skip=0;
$writit=0;
$num3= 3 * $numberofatoms;
$Matoms=0;
while (<FRAGJOB>) {

if(/Derivative of Dipole Moment Matrix/) {
$skip= $.+ 2;
}
if(/New transposed dipole drv/) {
$skip= $.+ 2;
}
if($skip == $.) {$writit=1};
if(/-----------------------------------------/) {$writit=0};
if( $writit == 1) {
$Matoms = $Matoms + 1;
}
if( $Matoms > $num3 ) {$writit=0};
if($writit == 1) {
 $line1="$_";
 chomp $line1;
 @line=split(/\s+/,$line1);
 $numline=$#line;
 print TMP "$line[$numline-2]        $line[$numline-1]        $line[$numline]\n";
}
}
close TMP;
close FRAGjob;
}

sub getenergy_qch {
 $TotEn=0;
 $secen=0;
 open(FRAGJOB,"<@_");
 while (<FRAGJOB>) {
  if (/Charge-charge energy/) {
   my @secline = split (/\s+/, $_);
   $secen = "$secline[4]";
  }
   if (/Total energy in the final/) {
    my @totline = split (/\s+/, $_);
    $ent = "$totline[9]";
    $TotEn = $ent - $secen;
   }
   if(/correlation energy/) {
     if(/Total/) {
      my @totline = split (/\s+/, $_);
      $entc = "$totline[6]";
      $TotEn = $TotEn + $entc;
     }
     if(/CCSD/) {
      my @totline = split (/\s+/, $_);
      my $last = $#totline;
      $entc = "$totline[$last]";
      $TotEn = $TotEn + $entc;
     }
    }
 }
}

sub getforces_qch   {
my $filename=$_[0];
open (FILE,"<$filename");
open (TMP,">tempfile$n");
 my $skip=0;
 my $writit=0;
 my $skip2=0;
 my $writit2=0;
 $numberofatoms=0;
 $firsttime=0;
 $firstgrad=0;
while (<FILE>) {
 if (/Full Analytical Gradient/) {
      $skip=$. + 1;
      $firstgrad=$firstgrad+1;
} else {
if (/Gradient of SCF/) {
      $skip=$. + 1;
      $firstgrad=$firstgrad+1;
}
}
 if ($. == $skip) {$writit=1};
 if (/Max gradient component/) {$writit=0};
 if (/Gradient time/) {$writit=0};

      if(( $writit == 1) && ($firstgrad == 1)) {
         print TMP "$_";
      }
 if (/\$molecule/) {
  $skip2=$.+1;
  $firsttime=$firsttime+1;
 }
 if ($. == $skip2) {$writit2=1};
 if (/\$end/) {$writit2=0};

      if( ($writit2 == 1) && ($firsttime == 1)) {
         $numberofatoms = $numberofatoms + 1;
      }
}
$numberofatoms=$numberofatoms-1;
close TMP;
close FILE;

# now parse the tempfile
open (TMP,">theseforces$n");
my $line1=`awk 'NR==1,NR==1 {print \$0}' tempfile$n`;
my @ar = split(/\s+/,"$line1");
my $nvar=$#ar;
my $i0=$numberofatoms / $nvar;
my @j0=split(/\./,"$i0");
my $i1=$j0[0];
my $mostatoms=$i1 * $nvar;
my $i2=$numberofatoms % $nvar;
for (my $j1=1; $j1 <= $i1; $j1++) {
 my $j3 = ($j1 - 1) * $nvar;
for (my $j2=2; $j2 <= 4; $j2++) {
 my $newline= ($j1 - 1) * 4 + $j2;
 my $line1=`awk 'NR=='$newline',NR=='$newline' {print \$0}' tempfile$n`;
 chomp $line1;
 my @grad=split(/\s+/,"$line1");
 for ($k1=2; $k1 <= $nvar+1; $k1++) {
  my $j4=$j3 + $k1 - 1;
  if ($j2 == 2) {$xx[$j4]=$grad[$k1]};
  if ($j2 == 3) {$yy[$j4]=$grad[$k1]};
  if ($j2 == 4) {$zz[$j4]=$grad[$k1]};
 }
}
}
$j1=$i1+1;
my $j3 = ($j1 - 1) * $nvar;
for (my $j2=2; $j2 <= 4; $j2++) {
 my $newline= ($j1 - 1) * 4 + $j2;
 my $line1=`awk 'NR=='$newline',NR=='$newline' {print \$0}' tempfile$n`;
 chomp $line1;
 my @grad=split(/\s+/,"$line1");
 for ($k1=2; $k1 <= $i2+1; $k1++) {
  my $j4=$j3 + $k1 - 1;
  if ($j2 == 2) {$xx[$j4]=$grad[$k1]};
  if ($j2 == 3) {$yy[$j4]=$grad[$k1]};
  if ($j2 == 4) {$zz[$j4]=$grad[$k1]};
 }
}

for ($j1=1; $j1 <= $numberofatoms; $j1++) {
$xx[$j1]=0-$xx[$j1];
$yy[$j1]=0-$yy[$j1];
$zz[$j1]=0-$zz[$j1];
 print TMP "$j1","       0     ","$xx[$j1]","       ","$yy[$j1]","       ","$zz[$j1]\n";
}
system("rm tempfile$n");

}

sub gethessian_qch {

my $filename="$_[0]";
# filename is an fchk file
my $numberofatoms=`awk '/Number of atoms/ {print \$5}' $filename`;
open (FH,"<$filename");
open (TMP,">hessian$n");
my $writit=0;
my $skip=0;
my $var=0;
while (<FH>) {
 if (/Cartesian Force Constants/) {
  $skip=$. + 1;
 }
  if(/Density/) {
    $writit=0;
  }
 if ($. == $skip) {$writit=1};
 if( $writit == 1) {
my $row="$_";
chomp $row;
my @hes=split(/\s+/,"$row");
my $lgh=$#hes;
for (my $ic=1;$ic <= $lgh; $ic++ ) {
$var=$var+1;
$hessian[$var]=$hes[$ic];
}
}
}

for (my $ic=1; $ic <= $var; $ic++) {
print TMP "$hessian[$ic]\n";
}
close FH;
close TMP;

}






