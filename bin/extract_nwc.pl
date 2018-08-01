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
  &getcoords_nwc("FRAG$ic.log");
  system("cat thesecoords$n >> FragDerivatives$n");
  &getenergy_nwc("FRAG$ic.log");
  open(FR,">>FragDerivatives$n");
  print FR "$TotEn\n";
  close FR;
  if ($deriv >= 1) {
   &getforces_nwc("FRAG$ic.log");
   system("cat theseforces$n >> FragDerivatives$n");
  }
  if ($deriv >= 2) {
   &gethessian_nwc("FRAG$ic",$natomf);
   system("cat hessian$n >> FragDerivatives$n");
   &getdipdr_nwc("FRAG$ic.log");
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
 &getenergy_nwc("ab.$ic1.$ic2.log");
 $TotE=$TotEn;
 $absign=`awk '/Isg_coeff/ {print \$2}' ab.$ic1.$ic2.com`;
 chomp $absign;
 if ($deriv >= 1) {
  open(TMP1,">>abforces$n");
  print TMP1 "$ic1","  ","$ic2\n";
  print TMP1 "$absign\n";
  &getforces_nwc("ab.$ic1.$ic2.log");
  $nat1=$numberofatoms;
  print TMP1 "$nat1\n";
  close TMP1;
  system("cat theseforces$n >> abforces$n");
 }
 if ($deriv >= 2) {
  open(TMP2,">>abhessians$n");
  print TMP2 "$ic1","  ","$ic2\n";
  print TMP2 "$absign\n";
  &gethessian_nwc("ab.$ic1.$ic2",$nat1);
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
 &getenergy_nwc("nb.$ic1.0.log");
 $TotE=$TotE - $TotEn;
 if ($deriv >= 1) {
  &getforces_nwc("nb.$ic1.0.log");
  $nat2=$numberofatoms;
  open(TMP1,">>abforces$n");
  print TMP1 "$ic1\n";
  print TMP1 "$nat2\n";
  close TMP1;
  system("cat theseforces$n >> abforces$n");
 }
 if ($deriv >= 2) {
  &gethessian_nwc("nb.$ic1.0",$nat2);
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
&getenergy_nwc("nb.$ic2.0.log");
 $TotE=$TotE - $TotEn;
 if ($deriv >= 1) {
  &getforces_nwc("nb.$ic2.0.log");
  $nat3=$numberofatoms;
  open(TMP1,">>abforces$n");
  print TMP1 "$ic2\n";
  print TMP1 "$nat3\n";
  close TMP1;
  system("cat theseforces$n >> abforces$n");
 }
 if ($deriv >= 2) {
  &gethessian_nwc("nb.$ic2.0",$nat3);
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

sub getcoords_nwc {
$filename="$_[0]";
 open(FRAGJOB,"<$filename");
 open(TMP,">thesecoords$n");
 my $skip=0;
 my $writit=0;
  while (<FRAGJOB>) {
if (/No.       Tag          Charge          X              Y              Z/) {
    $skip=$. + 2;
   }
 if ($. == $skip) {$writit=1};
 if(/Atomic Mass/)  {$writit=0};
 if( $writit == 1) {
 $line1="$_";
 @line=split(/\s+/,$line1);
 $b1=$#line;
 $b2=$b1-1;
 $b3=$b1-2;
 $b4=$b1-3;
 if($b1 >= 4) {
 $numele=int($line[$b4]);
 print TMP "0    $numele    0","   $line[$b3]","   $line[$b2]","   $line[$b1]\n";
}
  }
}
 close TMP;
 close FRAGJOB;
}

sub getdipdr_nwc {
# using the log file
open(FRAGJOB,"<$_[0]");
open(TMP,">thesedipdr$n");
$skipx=0;
$writitx=0;
$skipy=0;
$writity=0;
$skipz=0;
$writitz=0;
$num3=3 * $numberofatoms;
$nx=0;
$ny=0;
$nz=0;
while (<FRAGJOB>) {

if(/X vector of derivative dipole/) {
$skipx= $.+ 1;
}
if($skipx == $.) {$writitx=1};
if ($nx == $num3) {$writitx=0};
if($writitx == 1) {
$nx=$nx+1;
 $line1="$_";
 chomp $line1;
 @line=split(/\s+/,$line1);
 $numline=$#line;
 $numprint=$numline-2;
 print TMP "$line[$numprint]\n";
}
if(/Y vector of derivative dipole/) {
$skipy= $.+ 1;
}
if($skipy == $.) {$writity=1};
if ($ny == $num3) {$writity=0};
if($writity == 1) {
$ny=$ny+1;
 $line1="$_";
 chomp $line1;
 @line=split(/\s+/,$line1);
 $numline=$#line;
 $numprint=$numline-2;
 print TMP "$line[$numprint]\n";
}
if(/Z vector of derivative dipole/) {
$skipz= $.+ 1;
}
if($skipz == $.) {$writitz=1};
if ($nz == $num3) {$writitz=0};
if($writitz == 1) {
$nz=$nz+1;
 $line1="$_";
 chomp $line1;
 @line=split(/\s+/,$line1);
 $numline=$#line;
 $numprint=$numline-2;
 print TMP "$line[$numprint]\n";
}
}
close TMP;
close FRAGJOB;
}

sub getenergy_nwc  {
$filename="$_[0]";
$tce = `grep TCE $filename | wc -l`;
chomp($tce);
$tce=~ s/^\s+|\s+$//g;
if ($tce == 0) {
 $scf=`grep 'Total SCF energy' $filename | wc -l`;
 chomp($scf);
 $scf=~ s/^\s+|\s+$//g;
 if( $scf > 0 ) {
 system("grep 'Total SCF energy' $filename > ejunk$n");
 }
 $dft=`grep 'Total DFT energy' $filename | wc -l`;
 chomp($dft);
 $dft=~ s/^\s+|\s+$//g;
 if( $dft > 0 ) {
 system("grep 'Total DFT energy' $filename > ejunk$n");
 }
 $line1=`awk 'NR==1,NR==1 {print \$0}' ejunk$n`;
    @line=split(/\s+/,$line1);
    $size=$#line;
    $TotEn=$line[$size];
 $ans=`grep Correlation $filename | wc -l`;
 chomp($ans);
 $ans=~ s/^\s+|\s+$//g;
 if ( $ans > 0 ) {
  $corr=`awk '/Correlation energy/ {print \$3}' $filename`;
  chomp($corr);
  $TotEn=$TotEn + $corr;
 }
}else{
 $line1=`grep total $filename | grep energy | head -1`;
 @line=split(/\s+/,$line1);
 $size=$#line;
 $TotEn=$line[$size];
 }
$first=0;
open(FRAGJOB,"<$filename");
while (<FRAGJOB>) {
 if(/CCSD corr. energy/) {
  if( $first == 0 ) {
    my @totline = split (/\s+/, $_);
    my $last = $#totline;
    $entc = "$totline[$last]";
    $TotEn = $TotEn + $entc;
    $first=1;
  }
 }
 if(/\(T\) corr. energy/) {
    my @totline = split (/\s+/, $_);
    my $last = $#totline;
    $entc = "$totline[$last]";
    $TotEn = $TotEn + $entc;
 }
}
close FRAGJOB;
system("rm -f ejunk$n");
}

sub getforces_nwc {
$numberofatoms=0;
 $skip=0;
 $writit=0;
 open(TMP,">theseforces$n");
 open(FRAGJOB,"<@_");
 while (<FRAGJOB>) {
if (/atom               coordinates                        gradient/) {
 $skip=$. + 2;
 }
if ($. == $skip) {$writit=1};

if( $writit == 1) {
$thisline="$_";
chomp $thisline;
my @line=split(/\s+/,"$thisline");
my $val=$#line;
if ($val < 8) {
 close TMP;
 close FRAGJOB;
 return;
}
$numberofatoms=$numberofatoms+1;
$line[6]=0 - $line[6];
$line[7]=0 - $line[7];
$line[8]=0 - $line[8];
print TMP " 0    ","0    ","$line[6]   ","$line[7]   ","$line[8]\n";
}
}
 close TMP;
 close FRAGJOB;
}

sub gethessian_nwc  {
$file1="$_[0]";
system("cp $file1.hess hessian$n");
#@file=split(/\./,"$file1");
##system("cp $file[0].hess hessian");
}






