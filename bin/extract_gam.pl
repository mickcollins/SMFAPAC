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
  &getcoords_gam("FRAG$ic.log");
  system("cat thesecoords$n >> FragDerivatives$n");
  &getenergy_gam("FRAG$ic.log");
  open(FR,">>FragDerivatives$n");
  print FR "$TotEn\n";
  close FR;
  if ($deriv >= 1) {
   &getforces_gam("FRAG$ic.dat");
   system("cat theseforces$n >> FragDerivatives$n");
  }
  if ($deriv >= 2) {
   &gethessian_gam("FRAG$ic.log",$natomf);
   system("cat hessian$n >> FragDerivatives$n");
   &getdipdr_gam("FRAG$ic.log");
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
 &getenergy_gam("ab.$ic1.$ic2.log");
 $TotE=$TotEn;
 $absign=`awk '/Isg_coeff/ {print \$2}' ab.$ic1.$ic2.com`;
 chomp $absign;
 if ($deriv >= 1) {
  open(TMP1,">>abforces$n");
  print TMP1 "$ic1","  ","$ic2\n";
  print TMP1 "$absign\n";
  &getforces_gam("ab.$ic1.$ic2.dat");
  $nat1=$numberofatoms;
  print TMP1 "$nat1\n";
  close TMP1;
  system("cat theseforces$n >> abforces$n");
 }
 if ($deriv >= 2) {
  open(TMP2,">>abhessians$n");
  print TMP2 "$ic1","  ","$ic2\n";
  print TMP2 "$absign\n";
  &gethessian_gam("ab.$ic1.$ic2.log",$nat1);
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
 &getenergy_gam("nb.$ic1.0.log");
 $TotE=$TotE - $TotEn;
 if ($deriv >= 1) {
  &getforces_gam("nb.$ic1.0.dat");
  $nat2=$numberofatoms;
  open(TMP1,">>abforces$n");
  print TMP1 "$ic1\n";
  print TMP1 "$nat2\n";
  close TMP1;
  system("cat theseforces$n >> abforces$n");
 }
 if ($deriv >= 2) {
  &gethessian_gam("nb.$ic1.0.log",$nat2);
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
&getenergy_gam("nb.$ic2.0.log");
 $TotE=$TotE - $TotEn;
 if ($deriv >= 1) {
  &getforces_gam("nb.$ic2.0.dat");
  $nat3=$numberofatoms;
  open(TMP1,">>abforces$n");
  print TMP1 "$ic2\n";
  print TMP1 "$nat3\n";
  close TMP1;
  system("cat theseforces$n >> abforces$n");
 }
 if ($deriv >= 2) {
  &gethessian_gam("nb.$ic2.0.log",$nat3);
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


sub getcoords_gam {

 open(FRAGJOB,"<$_[0]");
 open(TMP,">thesecoords$n");

$skip=0;
$writit=0;
while (<FRAGJOB>) {

if(/CHARGE         X                   Y                   Z/) {
$skip= $.+ 1;
}
if($skip == $.) {$writit=1};
if(/INTERNUCLEAR DISTANCES/) {$writit=0};
if($writit == 1) {
if($_ eq "\n") {
}else{
 $line1="$_";
 chomp $line1;
 @line=split(/\s+/,$line1);
 $numele=int($line[2]);
 $line[2]=$line[3] / 1.8897259886;
 $line[3]=$line[4] / 1.8897259886;
 $line[4]=$line[5] / 1.8897259886;

 print TMP "0    $numele    0","   $line[2]","   $line[3]","   $line[4]\n";
}
}
}
close TMP;
close FRAGJOB;
}

sub getdipdr_gam {
open(FRAGJOB,"<$_[0]");
open(TMP,">thesedipdr$n");
$skip=0;
$writit=0;
while (<FRAGJOB>) {

if(/ATOM                 MU-X           MU-Y           MU-Z/) {
$skip= $.+ 1;
}
if($skip == $.) {$writit=1};
if($_ eq "\n") {$writit=0};
if($writit == 1) {
 $line1="$_";
 chomp $line1;
 @line=split(/\s+/,$line1);
 $numline=$#line;
 my $dx=$line[$numline-2] / 4.80289;
 my $dy=$line[$numline-1] / 4.80289;
 my $dz=$line[$numline] / 4.80289;
 print TMP "$dx      $dy      $dz\n";
}
}
close TMP;
close FRAGjob;
}

sub getenergy_gam {
$filename="$_[0]";
open(TMP,"<$filename");
while (<TMP>) {
if(/TOTAL ENERGY =/) {
 $line1="$_";
 chomp $line1;
 @line=split(/\s+/,$line1);
 my $ic=$#line;
 $TotEn=$line[$ic];
}
if(/E\(MP/) {
 $line1="$_";
 chomp $line1;
 @line=split(/\s+/,$line1);
 my $ic=$#line;
 $TotEn=$line[$ic];
}
if(/E\(CC/) {
 $line1="$_";
 chomp $line1;
 @line=split(/\s+/,$line1);
 my $ic=$#line;
 $TotEn=$line[$ic];
}
if(/E\(   CC/) {
 $line1="$_";
 chomp $line1;
 @line=split(/\s+/,$line1);
 my $ic=$#line;
 $TotEn=$line[$ic];
}
}

}

sub getforces_gam {
$filename="$_[0]";
$numberofatoms=0;
 $skip=0;
 $writit=0;
 open(TMP,">theseforces$n");
 open(FRAGJOB,"<$filename");
 while (<FRAGJOB>) {
 if (/GRAD /) {
 $skip=$. + 2;
 }
 if ($skip == $.) {$writit=1};
 if(/END/) {$writit=0};
 if($writit == 1) {
 $thisline="$_";
 chomp $thisline;
 my @line=split(/\s+/,"$thisline");
 $line[2]=0 - $line[2];
 $line[3]=0 - $line[3];
 $line[4]=0 - $line[4];
print TMP " 0    ","0    ","$line[2]   ","$line[3]   ","$line[4]\n";
$numberofatoms=$numberofatoms+1;
 }
 }
close TMP;
close FRAGJOB;
}

sub gethessian_gam {
  open(FRAGJOB,"<$_[0]");
 $numberofatoms=$_[1];
 $skip=0;
 $writit=0;
 open(TMP,">tfile$n");
 print TMP "$numberofatoms\n";
 while (<FRAGJOB>) {
 if(/X        Y        Z/) {
  $skip=$. + 1;
 }
 if($. == $skip) {$writit=1};
 if($_ eq "\n") {$writit=0};
 if($writit == 1) {
 $line1="$_";
 chomp $line1;
print TMP "$line1\n";
}


}
close TMP;
close FRAGJOB;

system("$EXEDIR/convertgamhessian < tfile$n > hessian$n");

}





