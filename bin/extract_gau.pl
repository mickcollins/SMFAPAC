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
  &getcoords_gau("FRAG$ic.log");
  system("cat thesecoords$n >> FragDerivatives$n");
  &getenergy_gau("FRAG$ic.log");
  open(FR,">>FragDerivatives$n");
  print FR "$TotEn\n";
  close FR;
  if ($deriv >= 1) {
   &getforces_gau("FRAG$ic.log");
   system("cat theseforces$n >> FragDerivatives$n");
  }
  if ($deriv >= 2) {
   &gethessian_gau("FRAG$ic.log",$natomf);
   system("cat hessian$n >> FragDerivatives$n");
   &getdipdr_gau("FRAG$ic.fchk");
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
 &getenergy_gau("ab.$ic1.$ic2.log");
 $TotE=$TotEn;
 $absign=`awk '/Isg_coeff/ {print \$2}' ab.$ic1.$ic2.com`;
 chomp $absign;
 if ($deriv >= 1) {
  open(TMP1,">>abforces$n");
  print TMP1 "$ic1","  ","$ic2\n";
  print TMP1 "$absign\n";
  &getforces_gau("ab.$ic1.$ic2.log");
  $nat1=$numberofatoms;
  print TMP1 "$nat1\n";
  close TMP1;
  system("cat theseforces$n >> abforces$n");
 }
 if ($deriv >= 2) {
  open(TMP2,">>abhessians$n");
  print TMP2 "$ic1","  ","$ic2\n";
  print TMP2 "$absign\n";
  &gethessian_gau("ab.$ic1.$ic2.log",$nat1);
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
 &getenergy_gau("nb.$ic1.0.log");
 $TotE=$TotE - $TotEn;
 if ($deriv >= 1) {
  &getforces_gau("nb.$ic1.0.log");
  $nat2=$numberofatoms;
  open(TMP1,">>abforces$n");
  print TMP1 "$ic1\n";
  print TMP1 "$nat2\n";
  close TMP1;
  system("cat theseforces$n >> abforces$n");
 }
 if ($deriv >= 2) {
  &gethessian_gau("nb.$ic1.0.log",$nat2);
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
&getenergy_gau("nb.$ic2.0.log");
 $TotE=$TotE - $TotEn;
 if ($deriv >= 1) {
  &getforces_gau("nb.$ic2.0.log");
  $nat3=$numberofatoms;
  open(TMP1,">>abforces$n");
  print TMP1 "$ic2\n";
  print TMP1 "$nat3\n";
  close TMP1;
  system("cat theseforces$n >> abforces$n");
 }
 if ($deriv >= 2) {
  &gethessian_gau("nb.$ic2.0.log",$nat3);
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





sub getcoords_gau {
 open(FRAGJOB,"<$_[0]");
 open(TMP,">thesecoords$n");
 my $skip=0;
 my $writit=0;
  while (<FRAGJOB>) {
   if (/Input orientation/) {
    $skip=$. + 5;
   }  else {
    if (/Z-Matrix orientation/) {
     $skip=$. + 5;
    }
   }
 if ($. == $skip) {$writit=1};
 if (/-----------------/) {$writit=0};

      if( $writit == 1) {
         print TMP "$_";
      }
  }
 close TMP;
 close FRAGJOB;
}

sub getdipdr_gau {
# the input file is a FRAG.*.fchk file
open(FRAGJOB,"<$_[0]");
open(TMP,">thesedipdr$n");
$skip=0;
$writit=0;
while (<FRAGJOB>) {

if(/Dipole Derivatives/) {
$skip= $.+ 1;
}
if($skip == $.) {$writit=1};
if(/Polarizability/) {$writit=0};
if($writit == 1) {
 $line1="$_";
 chomp $line1;
 print TMP "$line1\n";
}
}
close TMP;
close FRAGjob;
}

sub getenergy_gau {
 $TotEn=0;
 $firstE2=1;
 $firstSCF=1;
 open(FRAGJOB,"<@_");
 while (<FRAGJOB>) {
  if (/Self energy of the charges/) {
   my @secline = split (/\s+/, $_);
   my $secen = "$secline[7]";
   $TotEn = $TotEn - $secen;
  }
  if(/Convergence criterion not met/) {
   $firstSCF=0;
  }
  if(/convergence achieved/) {
   $firstSCF=1;
  }
  if (/SCF Done:/)  {
   my @scfline = split (/\s+/, $_);
   my $scfen = "$scfline[5]";
   $TotEn = $TotEn + $scfen*$firstSCF;
  }
  if (/E2 =/)  {
   my @mp2line = split (/\s+/, $_);
   my $mp2en = "$mp2line[3]";
   $mp2en =~ s/[dD]/E/g;
   $TotEn = $TotEn + $mp2en*$firstE2;
   $firstE2=0;
  }
  if(/CCSD\(T\)=/) {
   if( $firstCC == 1 ) {
    my @line = split (/\s+/, $_);
    $TotEn = $line[2];
    $firstCC = 0;
   }
  }
  if(/E\(CORR\)/) {
   $ccsd="$_";
  }
 }

  close FRAGJOB;
 if ( $ccsd ne "nothing" ) {
  if($firstCC == 1) {
  my @line = split (/\s+/,$ccsd);
  $TotEn = "$line[4]";
  }
 }

 open(TMP,">thisenergy$n");
 print TMP "$TotEn\n";
 close TMP;
}

sub getforces_gau {

 $skip=0;
 $writit=0;
 $numberofatoms=0;
 open(TMP,">theseforces$n");
 open(FRAGJOB,"<@_");
 while (<FRAGJOB>) {
  if (/Forces \(Hartrees\/Bohr\)/) {
   $skip=$. + 3;
  }
if ($. == $skip) {$writit=1};
if (/-----------------/) {$writit=0};

      if( $writit == 1) {
         print TMP "$_";
         $numberofatoms = $numberofatoms + 1;
      }
 }
 close TMP;
 close FRAGJOB;
}

sub gethessian_gau {

 $skip=0;
 $writit=0;
 open(TMP,">archive$n");
 open(FRAGJOB,"<$_[0]");
 my $natomf=$_[1];
 while (<FRAGJOB>)  {
  if(/GINC/) {
 $skip=$.;
 }
 if ($. == $skip) {$writit=1};

      if( $writit == 1) {
         print TMP "$_";
      }

 if (/@/) {$writit=0};
 }
 close TMP;
system("$CODEDIR/striparchive_gau_par archive$n $n");
my $n3 = 3 * $natomf;
my $n4 = $n3 * $n3 + $n3;
my $nsec = $n4 / 2;
my $ntail = $nsec + $n3 + 1;
system("tail -$ntail rub2_$n > rub1_$n");
system("head -$nsec rub1_$n > hessian$n");
 close FRAGJOB;
}




