sub readgauparams {

if ( -s "xyzFILENAME" ) {
$file=`awk 'NR==1,NR==1 {print \$1}' xyzFILENAME`;
}else{
print "You must first enter the filename for the Cartesian coordinates";
system ("sleep 5");
return 1;
}

system("$EXEDIR/distinctlabels < $file ");
system("mv labels labels0");

print "                        Variables need by SMFA and GAUSSIAN\n";
print "\n";

$deriv=0;
print "Enter the type of calculation as an integer:\n";
print "Energy    (0)\n"; 
print "Force     (1)\n";
print "Frequency (2)\n";
print "Optimize  (3)\n";
print "Find TS   (4)\n";
print "Scan      (5)\n";
$deriv1=<STDIN>;
chomp $deriv1;
if($deriv1 == 1) {$deriv=1};
if($deriv1 == 2) {$deriv=2};
if($deriv1 == 3) {$deriv=3};
if($deriv1 == 4) {$deriv=4};
if($deriv1 == 5) {$deriv=5};
if($deriv==0){$prop=" "};
if($deriv==1){$prop="FORCE"};
if($deriv==2){$prop="FREQ"};
if($deriv==3){$prop="OPT"};
if($deriv==4){$prop="OPT"};
if($deriv==5){$prop="OPT"};

print "\n";

print "Enter the calculation method [eg HF or MP2 or other]\n";
$method=<STDIN>;
chomp $method;
$method=uc($method);
print "\n";
print "Enter the %mem value, eg %mem=500mb, or hit RETURN for the GAUSSIAN default value on your system\n";
$mem=<STDIN>;
if ($mem eq "\n") {
}else{
chomp $mem;
}
print "\n";
print "Does this method account for electronic correlation leading to dispersion?\n";
print "(See the user's manual). Do you want to account for dispersion at long range?\n";
print "Answer Y or N\n";
$disp=<STDIN>;
chomp $disp;
$disp=uc($disp);
if($disp ne "Y") {
$disp="N";
}
print "\n";
print "GAUSSIAN allows a different basis set for each chemical element.\n";
print "Do you want all elements to have the same basis (Y or N) ? \n";
$ans=<STDIN>;
chomp $ans;
$ans=uc($ans);
if($ans ne "N") {$ans="Y"};
if ($ans eq "Y") {
print "Enter the basis for all atoms\n";
$basis=<STDIN>;
chomp $basis;
system("rm -f GAUbasis");
}else{
print "Enter the basis (check availability in the GAUSSIAN manual)\n";
print "for the following elements\n";
open (TMP,"<labels0");
system("rm -f GAUbasis");
open (TMP1,">GAUbasis");
while (<TMP>) {
my $lab="$_";
chomp $lab;
print "$lab\n";
my $abasis=<STDIN>;
chomp $abasis;
print TMP1 "$lab  0\n";
print TMP1 "$abasis\n";
print TMP1 "****\n";
}
close TMP;
close TMP1;
}

print "\n";
print "Enter any other keywords (optional), or hit RETURN:\n ";
print "Enter all keywords on one line (with appropriate spaces)\n";
$opts=<STDIN>;
chomp $opts;
$opts=uc($opts);
#print "\n";
#print "Enter the number of processors for each calculation:\n ";
#$nproc0=<STDIN>;
#chomp $nproc0;

&solventcharges;

if($deriv==3) {&optinput};
if($deriv==4) {&optinput};
if($deriv==4) {&TSinput};
if($deriv==5) {&scaninput};


      open TMP, ">IN_JOBTYPE";
      print TMP "Enter 0, 1, 2, 3, 4, 5 for energy, gradient, hessian, opt, TS or scan\n";
      print TMP "$deriv\n";
      close TMP;

      open TMP, ">ABSETTINGS";
print TMP "The quantum chemistry program package is\n";
       print TMP "$package\n";
print TMP "The electronic structure method is\n";
       print TMP "$method\n";
if ($ans eq "Y") {
 print TMP "The basis set is\n";
 print TMP "$basis\n";
}else{
print TMP "The basis gen sets are\n";
close TMP;
system("cat GAUbasis >> ABSETTINGS");
open TMP, ">>ABSETTINGS";
}
print TMP "The job type is energy (0), force (1), hessian (2), opt (3), opt=TS (4), scan (5)\n";
       print TMP "$deriv\n";
print TMP "Additional options are\n";
       print TMP "$opts\n";
print TMP "Is long range dispersion accounted for?\n";
       print TMP "$disp\n";
#print TMP "The number of processors for individual calculations is\n";
#       print TMP "$nproc0\n";
      close TMP;

if ($deriv >= 1 ) {
 $conv="SCF=Tight";
}else{
 $conv="";
}

$dens="density=current";
if ($method eq "CCSD\(T\)") {
 $dens="";
}
# Gaussian job header
       open G09IN, ">IN_G09";
#       # Gaussian job Model - method/basis
#             print G09IN "\%nproc=$nproc0\n";
print G09IN "%nproc=1\n";
if ($mem ne "\n") {
print G09IN "$mem\n";
}
if ($ans eq "Y") {
             print G09IN "\#p $method/$basis nosymm $dens $prop $conv $opts\n\n";
}else{
             print G09IN "\#p $method/gen nosymm $dens $prop  $conv $opts\n\n";
}
             close G09IN;
       open G09IN, ">IN_NPA";
$popnpa="pop=NPA";
# header for NPA charges
#             print G09IN "\%nproc=$nproc0\n";
print G09IN "%nproc=1\n";
if ($mem ne "\n") {
print G09IN "$mem\n";
}
if ($ans eq "Y") {
print G09IN "\#p $method/$basis charge nosymm $popnpa\n\n";
}else{
             print G09IN "\#p $method/gen charge nosymm $popnpa\n\n";
}
             close G09IN;

       open G09IN, ">IN_POLAR";
#             print G09IN "\%nproc=$nproc0\n";
print G09IN "%nproc=1\n";
if ($mem ne "\n") {
print G09IN "$mem\n";
}
if ($ans eq "Y") {
print G09IN "\#p $method/$basis polar nosymm\n\n";
}else{
print G09IN "\#p $method/gen polar nosymm\n\n";
}
             close G09IN;
# make a special IN_G09 for the initial geometry in an optimisation/ts/scan
       open G09IN, ">IN_G09_step0";
print G09IN "%nproc=1\n";
if ($mem ne "\n") {
print G09IN "$mem\n";
}
if ($ans eq "Y") {
             print G09IN "\#p HF/$basis FREQ nosymm density=current\n\n";
}else{
             print G09IN "\#p HF/gen FREQ nosymm density=current\n\n";
}
             close G09IN;


# make the IN_GDMA file
open(TMP,">IN_GDMA");
$den="SCF";
if($method eq "MP2") {$den="MP2"};
if($method eq "CCSD") {$den="CC"};
print TMP "Density $den\n";
print TMP "\n";
print TMP "Multipoles\n";
print TMP " Switch 0\n";
print TMP " Limit 4\n";
print TMP " Limit 4 H\n";
print TMP " Punch\n";
print TMP "Start\n";
print TMP "\n";
print TMP "Finish\n";
close TMP;
open(TMP,">IN_GDMA_step0");
$den="SCF";
print TMP "Density $den\n";
print TMP "\n";
print TMP "Multipoles\n";
print TMP " Switch 0\n";
print TMP " Limit 4\n";
print TMP " Limit 4 H\n";
print TMP " Punch\n";
print TMP "Start\n";
print TMP "\n";
print TMP "Finish\n";
close TMP;



}


sub extract_gau {

# get the energy, force and hessian from all re;evant ab initio jobs
# run by GAUSSIAN

system("clear");
system("stty echo");

# first the FRAG jobs
 my $nfrgs = `ls COORD* | wc -l`;
 open(FR,">FragDerivatives");
 print FR "$deriv\n";
 print FR "$nfrgs";
 close FR;
if($deriv >= 2) {
open(DR,">FragDipderivs");
 print DR "$nfrgs";
 close DR;
}

 for (my $ic=1; $ic<=$nfrgs; $ic++) {

my $mf=`wc -l COORD$ic`;
my @fields = split / /, $mf;
$natomf=$fields[0];
$natomf = $natomf - 1;
open(FR,">>FragDerivatives");
print FR "$natomf\n";
close FR;


&getcoords_gau("FRAG$ic.log");
system("cat thesecoords >> FragDerivatives");

&getenergy_gau("FRAG$ic.log");
open(FR,">>FragDerivatives");
print FR "$TotEn\n";
close FR;
  if ($deriv >= 1) {
   &getforces_gau("FRAG$ic.log");
   system("cat theseforces >> FragDerivatives");
  }
  if ($deriv >= 2) {
   &gethessian_gau("FRAG$ic.log",$natomf);
   system("cat hessian >> FragDerivatives");
   &getdipdr_gau("FRAG$ic.fchk");
   system("echo '$numberofatoms' >> FragDipderivs");
   system("cat thesedipdr >> FragDipderivs");
  }
 }
# now the ab initio nonbonded jobs
system("ls -1 ab.*.log > abfiles");
my $ab=`wc -l abfiles`;
my @fields = split / /, $ab;
my $nfragsab=$fields[0];
if ($deriv >= 1) {
 system("rm -f abforces");
 open(TMP1,">>abforces");
 print TMP1 "$nfragsab\n";
}
if ($deriv >= 2) {
 system("rm -f abhessians");
 open(TMP2,">>abhessians");
 print TMP2 "$nfragsab\n";
}

$filename='abfiles';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
my $AllTot = 0;
while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 my $ic2=$fields[2];
 &getenergy_gau("ab.$ic1.$ic2.log");
 my $TotE=$TotEn;
 &getenergy_gau("nb.$ic1.0.log");
 $TotE=$TotE - $TotEn;
 &getenergy_gau("nb.$ic2.0.log");
 $TotE=$TotE - $TotEn;
 my $absign=`awk '/Isg_coeff/ {print \$2}' ab.$ic1.$ic2.com`;
 chomp $absign;
 $AllTot = $AllTot + $absign * $TotE;

 if ($deriv >= 1) {
  print TMP1 "$ic1","  ","$ic2\n";
  print TMP1 "$absign\n";
  &getforces_gau("ab.$ic1.$ic2.log");
  print TMP1 "$numberofatoms\n";
  my $nat1=$numberofatoms;
  close TMP1;
  system("cat theseforces >> abforces");
  open(TMP1,">>abforces");
  print TMP1 "$ic1\n";
  &getforces_gau("nb.$ic1.0.log");
  print TMP1 "$numberofatoms\n";
  my $nat2=$numberofatoms;
  close TMP1;
  system("cat theseforces >> abforces");
  open(TMP1,">>abforces");
  print TMP1 "$ic2\n";
  &getforces_gau("nb.$ic2.0.log");
    print TMP1 "$numberofatoms\n";
    my $nat3=$numberofatoms;
  close TMP1;
  system("cat theseforces >> abforces");
  open(TMP1,">>abforces");


 if ($deriv >= 2) {
  print TMP2 "$ic1","  ","$ic2\n";
  print TMP2 "$absign\n";
  &gethessian_gau("ab.$ic1.$ic2.log",$nat1);
  print TMP2 "$nat1\n";
  close TMP2;
  system("cat hessian >> abhessians");
  open(TMP2,">>abhessians");
  print TMP2 "$ic1\n";
  &gethessian_gau("nb.$ic1.0.log",$nat2);
  print TMP2 "$nat2\n";
  close TMP2;
  system("cat hessian >> abhessians");
  open(TMP2,">>abhessians");
  print TMP2 "$ic2\n";
  &gethessian_gau("nb.$ic2.0.log",$nat3);
    print TMP2 "$nat3\n";
  close TMP2;
  system("cat hessian >> abhessians");
  open(TMP2,">>abhessians");
 }
# end the deriv=1 if
 }
}
close TMP1;
close TMP2;
open (TMP1,">NearInt_energy");
print TMP1 "$AllTot\n";
close TMP1;

}

sub getcoords_gau {
 open(FRAGJOB,"<$_[0]");
 open(TMP,">thesecoords");
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
open(TMP,">thesedipdr");
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

 open(TMP,">thisenergy");
 print TMP "$TotEn\n";
 close TMP;
}

sub getforces_gau {

 $skip=0;
 $writit=0;
 $numberofatoms=0;
 open(TMP,">theseforces");
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
 open(TMP,">archive");
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
system("$CODEDIR/striparchive_gau archive");
my $n3 = 3 * $natomf;
my $n4 = $n3 * $n3 + $n3;
my $nsec = $n4 / 2;
my $ntail = $nsec + $n3 + 1;
system("tail -$ntail rub2 > rub1");
system("head -$nsec rub1 > hessian");
 close FRAGJOB;
}

sub makenbgdma_gau {
$step=1;
if ( -s "STEPNUMBER" ) {
$step=`awk 'NR==1,NR==1 {print \$1}' STEPNUMBER`;
chomp($step);
$step=~ s/^\s+|\s+$//g;
}
if ($step == 0) {
open(GDMA,"IN_GDMA_step0") or die "Unable to open IN_GDMA";
}else{
open(GDMA,"IN_GDMA") or die "Unable to open IN_GDMA";
}
@GDMALines = <GDMA>;

$pwd=`pwd`;
chomp $pwd;

@All_jobs= `ls -lah nb*.fchk --format=single-column`;
         $Number_of_jobs=@All_jobs;
         foreach $job (@All_jobs){
            chomp $job;
            $prefix=$job;
            $prefix =~ s/\.[^.]*$//;
#            `formchk $job $prefix.fchk`;
            $gdma_file="$prefix-gdma.inp";
             if(-s "$gdma_file")
             {
             open (my $gdma, "<$gdma_file") or die "Unable to open $!";
             @Names=<$gdma>;
             close $gdma;
             }
             open (my $gdma, ">$gdma_file") or die "Unable to open $!";
             &make_gdma($gdma);

            `$EXEDIR/gdma <$gdma_file>$prefix-gdma.out`;
            `$EXEDIR/Stonepunch2cart_Matt <$prefix.punch>$prefix.cart`;
#            `cart2chg.pl $prefix.cart`;
          }
}

sub ReadPolar_gau  {
system("rm Polarisabilities.out");
open(my $OUT, ">Polarisabilities.out");

my $nfrgs = `ls nb.*.0-polar.log | wc -l`;
chomp $nfrgs;
for (my $ic=1; $ic<=$nfrgs; $ic++) {
$job="nb.$ic.0-polar.log";
&getpolar_gau("$job");
$iso=($Polar[1] + $Polar[3] + $Polar[6])/3;
print $OUT "$ic\t$Polar[1]\t$Polar[2]\t$Polar[3]\t$Polar[4]\t$Polar[5]\t$Polar[6]\t$iso\n";

}
}


sub getpolar_gau {
 $infile="$_[0]";

 open(FRAGJOB,"<$infile");

while (<FRAGJOB>)  {
 if (/Exact polarizability/) {
  my $line="$_";
  my @ans=split(/\s+/,"$line");
  $Polar[1]=$ans[3];
  $Polar[2]=$ans[4];
  $Polar[3]=$ans[5];
  $Polar[4]=$ans[6];
  $Polar[5]=$ans[7];
  $Polar[6]=$ans[8];
 }
}
}

sub makehyperpolar_gau {
system("clear");
system("stty sane");

$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
$n++;
}
$n=$n-1;
#$nproc0=$store[$n];

open(TMP,"<IN_JOBTYPE");
while (<TMP>) {
if ($. == 2) {$deriv="$_"};
}
chomp $deriv;

$qmprog=$store[1];

if( ($qmprog != 2) || ($deriv != 2)) {
print "Hyperpolarizabilities are only available for GAUSSIAN following a frequency calculation\n";
system("sleep 5");
return;
}

open(TMP,">>OUT_SMFA");
print TMP "\n";
print TMP "The total molecular hyperpolarizability is estimated from a sum over\n";
print TMP "the hyperpolarizabilities of the fragments.\n";
print TMP "\n";

my $nfrgs = `ls COORD* | wc -l`;
chomp $nfrgs;
open(SG,"<signs.out");
while(<SG>) {
$sign[$.]=$_;
chomp $sign[$.];
$sign[$.]=ltrim($sign[$.]);
}
close SG;
for (my $ic=1;$ic <= 10;$ic++) {
$tothyper[$ic]=0;
}
for (my $ic=1;$ic <= $nfrgs;$ic++) {
&gethyperpolar_gau("FRAG$ic.fchk");
for (my $ie=1;$ie <= 10;$ie++) {
 $tothyper[$ie]=$tothyper[$ie]+$sign[$ic]*$HypPolar[$ie];
}
}
print TMP "      xxx        xyy        xyy        yyy        xxz\n";
printf TMP "%10.3f %10.3f %10.3f %10.3f %10.3f\n",$tothyper[1],$tothyper[2],$tothyper[3],$tothyper[4],$tothyper[5];
print TMP "\n";
print TMP "      xyz        yyz        xzz        yzz        zzz\n";
printf TMP "%10.3f %10.3f %10.3f %10.3f %10.3f\n",$tothyper[6],$tothyper[7],$tothyper[8],$tothyper[9],$tothyper[10];
print TMP "\n";
close TMP;
print "The hyperpolarizability tensor has been appended to OUT_SMFA\n";
system("sleep 4");

}

sub gethyperpolar_gau {
 $infile="$_[0]";
 $writit=0;
$linenum1=0;
$linenum2=0;
 open(FRAGJOB,"<$infile");
while (<FRAGJOB>)  {
if (/HyperPolarizability/)  {
  $writit=1;
  $linenum1=$.+1;
  $linenum2=$linenum1+1;
}
if($linenum1 == $.) {
  my $line="$_";
  my @ans=split(/\s+/,"$line");
  my $l5=$#ans;
  my $l4=$l5-1;
  my $l3=$l4-1;
  my $l2=$l3-1;
  my $l1=$l2-1;
  $HypPolar[1]=$ans[$l1];
  $HypPolar[2]=$ans[$l2];
  $HypPolar[3]=$ans[$l3];
  $HypPolar[4]=$ans[$l4];
  $HypPolar[5]=$ans[$l5];
}
if($linenum2 == $.) {
  my $line="$_";
  my @ans=split(/\s+/,"$line");
  my $l5=$#ans;
  my $l4=$l5-1;
  my $l3=$l4-1;
  my $l2=$l3-1;
  my $l1=$l2-1;
  $HypPolar[6]=$ans[$l1];
  $HypPolar[7]=$ans[$l2];
  $HypPolar[8]=$ans[$l3];
  $HypPolar[9]=$ans[$l4];
  $HypPolar[10]=$ans[$l5];
}

}
}


sub SMFArungau {

  my $nfrgs = `ls COORD* | wc -l`;
 chomp $nfrgs;
 for (my $ic=1; $ic<=$nfrgs; $ic++)
  {
  system("g09 < FRAG$ic.com > FRAG$ic.log");
 $ok=`grep 'Normal termination' FRAG$ic.log | wc -l`;
 chomp($ok);
 $ok=~ s/^\s+|\s+$//g;
 if ( $ok == 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation failed for FRAG$ic\n";
  close TMP;
  exit(0);
 }
 if($deriv >= 2) {
 system("formchk FRAG$ic.chk");
 }
  }

 system("ls -1 ab*.com > abfiles");
$filename='abfiles';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 my $ic2=$fields[2];
 system("g09 < ab.$ic1.$ic2.com > ab.$ic1.$ic2.log");
 $ok=`grep 'Normal termination' ab.$ic1.$ic2.log | wc -l`;
 chomp($ok);
 $ok=~ s/^\s+|\s+$//g;
 if ( $ok == 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation failed for ab.$ic1.$ic2\n";
  close TMP;
  exit(0);
 }
 }
close $fh;

 system("ls -1 nb*.0.com > nbfiles");
$filename='nbfiles';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 system("g09 < nb.$ic1.0.com > nb.$ic1.0.log");
 $ok=`grep 'Normal termination' nb.$ic1.0.log | wc -l`;
 chomp($ok);
 $ok=~ s/^\s+|\s+$//g;
 if ( $ok == 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation failed for nb.$ic1.0\n";
  close TMP;
  exit(0);
 }
 system("formchk nb.$ic1.0.chk nb.$ic1.0.fchk");
 }
close $fh;

$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
$n++;
}
$n=$n-1;
#$nproc0=$store[$n];
close $fh;

for (my $ic=0; $ic<=$n; $ic++) {
 if ("$store[$ic]" eq "Is long range dispersion accounted for?")  {
  $ie=$ic+1;
  $disp=$store[$ie];
 }
}

    $nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;
if ($nchgs eq 0) {

 system("ls -1 nb*.0-polar.com > nbfilesp");
$filename='nbfilesp';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 system("g09 < nb.$ic1.0-polar.com > nb.$ic1.0-polar.log");
 $ok=`grep 'Normal termination' nb.$ic1.0-polar.log | wc -l`;
 chomp($ok);
 $ok=~ s/^\s+|\s+$//g;
 if ( $ok == 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation failed for nb.$ic1.0-polar\n";
  close TMP;
  exit(0);
 }

 }
close $fh;
} else {
if ("$disp" eq "Y") {
 system("ls -1 nb*.0-polar.com > nbfilesp");
$filename='nbfilesp';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 system("g09 < nb.$ic1.0-polar.com > nb.$ic1.0-polar.log");
 $ok=`grep 'Normal termination' nb.$ic1.0-polar.log | wc -l`;
 chomp($ok);
 $ok=~ s/^\s+|\s+$//g;
 if ( $ok == 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation failed for nb.$ic1.0-polar\n";
  close TMP;
  exit(0);
 }

 }
close $fh;
}
}

}

sub Lev0_chg_MAC_iter_gau {

system("clear");

$iter=3;
$count=1;
while ($count <= $iter) {

&SMFAgaunpa_MAC;
&extractch_gau;

      $count++;
}

      return 1;
    }

sub extractch_gau {
      my @jobs = `ls charge.*.coord`;
      foreach my $chg_coord (@jobs) {
      @fields=split(/\./,$chg_coord);
      $ic=$fields[1];

      open (CHGFIL,"$chg_coord");
      $n=0;
      while (<CHGFIL>) {
      chomp;
      if ($. == 1) {
      ($dum,$chg,$mlt) = split /\s+/, $_;
      } else {
      ($symb[$n],$x[$n],$y[$n],$z[$n])=split(/\s+/,$_);
      $n++;
      }
      $Nat = $n;
      }
      close CHGFIL;


$start=0;
$end=0;

      open (NPA_IN,"charge.$ic.log");
      $n=0;
      while (<NPA_IN>) {
      if (/Summary of Natural Population Analysis/){
      $start=$. + 5;
      $end=$start+$Nat;
      }
      if (($. > $start) && ($. <= $end)) {
      ($dum,$dum,$dum,$chgarray[$n]) = split (/\s+/, $_);
      $n++;
      }
      }
      close NPA_IN;
      open (NPA_OUT,">charge.$ic.npa");
      print NPA_OUT "charge.$ic.coord\n";
      print NPA_OUT "$Nat\n";

      for ($n=0; $n < $Nat; $n++) {
      print NPA_OUT "$x[$n] $y[$n] $z[$n] $chgarray[$n]\n";
      }
      close NPA_OUT;

      system("$EXEDIR/adjustcaponly < charge.$ic.npa");

      }
# end sub
}

sub SMFAgaunpa_MAC {

&SMFAmkchinputs_gau;
&runchseq_gau;

}

sub SMFAmkchinputs_gau {

$nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;

 my @jobs = `ls charge.*.coord`;
 foreach my $chg_coord (@jobs) {
  @fields=split(/\./,$chg_coord);
  $ic=$fields[1];
  system("cp IN_NPA charge.$ic.com");

  system("echo 'npa job' >> charge.$ic.com");
  system("echo ' ' >> charge.$ic.com");
  system("cat charge.$ic.coord >> charge.$ic.com");
  system("echo ' ' >> charge.$ic.com");
if($nchgs > 1) {
  open (CHID,"chL0.$ic");
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $m1 = &ltrim($n1);
   system("cat charge.$m1.grp >> charge.$ic.com");
   }
  system("echo ' ' >> charge.$ic.com");
}
$ans=`grep gen IN_G09 | wc -l`;
chomp $ans;
  if ($ans ne 0) {
# get only the basis for the elements present
open (TMP,">>charge.$ic.com");
my $len1=`wc -l GAUbasis`;
chomp $len1;
my @len=split(/\s+/,"$len1");
my $length=$len[0];
for (my $ig = 1; $ig <= $length; $ig +=3) {
$ele=`awk 'NR=='$ig',NR=='$ig' {print \$1}' GAUbasis`;
chomp $ele;
$here=`grep "$ele" "charge.$ic.coord" | wc -l`;
chomp $here;
if ($here > 0) {
my $line1=`awk 'NR=='$ig',NR=='$ig' {print \$0}' GAUbasis`;
my $line2=`awk 'NR=='$ig+1',NR=='$ig+1' {print \$0}' GAUbasis`;
my $line3=`awk 'NR=='$ig+2',NR=='$ig+2' {print \$0}' GAUbasis`;
chomp $line1;
chomp $line2;
chomp $line3;
print TMP "$line1\n";
print TMP "$line2\n";
print TMP "$line3\n";
}
}
print TMP " \n";
close TMP;

  }
  }
# end sub
}

sub runchseq_gau {

$nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;

 my @jobs = `ls charge.*.coord`;
 foreach my $chg_coord (@jobs) {
  @fields=split(/\./,$chg_coord);
  $ic=$fields[1];
  system("g09 < charge.$ic.com >  charge.$ic.log");
 $ok=`grep 'Normal termination' charge.$ic.log | wc -l`;
 chomp($ok);
 $ok=~ s/^\s+|\s+$//g;
 if ( $ok == 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation failed for charge.$ic\n";
  close TMP;
  exit(0);
 }
}
# end sub
}


sub SMFAgauinputs_MAC  {

    $nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;

my $pwd=`pwd`;
chomp $pwd;

#open(IN,">INLIST");

# start with the main fragments
  my $nfrgs = `ls COORD* | wc -l`;
 chomp $nfrgs;
 for (my $ic=1; $ic<=$nfrgs; $ic++)
  {
  my $my_coord = "COORD${ic}";
# find out if there are any external charges
  $ffile[0]=0;
  if( -s "chFRAG.$ic") {
  $efile=`wc -l chFRAG.$ic`;
  chomp $efile;
  @ffile=split(/\s+/,$efile);
  }
# if necessary add charge to the G09 header
  system("cp IN_G09 IN_G09_temp");
  if($ffile[0] ne "0") {
   system("sed 's/nosymm/charge nosymm/' IN_G09_temp > IN_G09_current");
  }else{
   system("cp IN_G09 IN_G09_current");
  }
  if($deriv >= 2) {
  system("echo '%chk=$pwd/FRAG$ic.chk' > FRAG$ic.com");
  }
  system("cat IN_G09_current >> FRAG$ic.com");
  system("echo 'Bonded fragment', $ic, ' with any necessary charges\n' >> FRAG$ic.com");
  system("cat $my_coord >> FRAG$ic.com");
  system("echo ' ' >> FRAG$ic.com");
    if($nchgs >= 1){
  open (CHID,"chFRAG.$ic");
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $m1 = ltrim($n1);
   system("cat charge.$m1.grp >> FRAG$ic.com");
   }
  close CHID;
    }
  system("echo ' ' >> FRAG$ic.com");
# now add the gen basis if necessary
$ans=`grep gen IN_G09 | wc -l`;
chomp $ans;
  if ($ans ne 0) {
# get only the basis for the elements present
open (TMP,">>FRAG$ic.com");
my $len1=`wc -l GAUbasis`;
chomp $len1;
my @len=split(/\s+/,"$len1");
my $length=$len[0];
for (my $ig = 1; $ig <= $length; $ig +=3) {
$ele=`awk 'NR=='$ig',NR=='$ig' {print \$1}' GAUbasis`;
chomp $ele;
$here=`grep "$ele" "COORD$ic" | wc -l`;
chomp $here;
if ($here > 0) {
my $line1=`awk 'NR=='$ig',NR=='$ig' {print \$0}' GAUbasis`;
my $line2=`awk 'NR=='$ig+1',NR=='$ig+1' {print \$0}' GAUbasis`;
my $line3=`awk 'NR=='$ig+2',NR=='$ig+2' {print \$0}' GAUbasis`;
chomp $line1;
chomp $line2;
chomp $line3;
print TMP "$line1\n";
print TMP "$line2\n";
print TMP "$line3\n";
}
}
print TMP " \n";
close TMP;

  }
#print IN "FRAG$ic.com\n";
    }

# now the ab.*.* files
 system("ls -1 ab.*.com > abfiles");
$filename='abfiles';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 my $ic2=$fields[2];
# find out if there are any external charges
  $ffile[0]=0;
  if( -s "chab.$ic1.$ic2") {
  $efile=`wc -l chab.$ic1.$ic2`;
  chomp $efile;
  @ffile=split(/\s+/,$efile);
  }
# if necessary add charge to the G09 header
  system("cp IN_G09 IN_G09_temp");
  if($ffile[0] ne "0") {
   system("sed 's/nosymm/charge nosymm/' IN_G09_temp > IN_G09_current");
  }else{
   system("cp IN_G09 IN_G09_current");
  }

 system("cp IN_G09_current tmpfile");
 system("cat ab.$ic1.$ic2.com >> tmpfile");
# system("echo ' ' >> tmpfile");
    if($nchgs >= 1){
  open (CHID,"chab.$ic1.$ic2");
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $m1 = ltrim($n1);
   system("cat charge.$m1.grp >> tmpfile");
   }
  close CHID;
  system("echo ' ' >> tmpfile");
    }
# now add the gen basis if necessary
  if ($ans ne 0) {
# find the length of the ab file to avoid the first line
my $len1=`wc -l ab.$ic1.$ic2.com`;
chomp $len1;
my @len=split(/\s+/,"$len1");
my $lengthab=$len[0];
$lengthab=$lengthab-1;
system("tail -$lengthab ab.$ic1.$ic2.com > tfile");
open (TMP,">>tmpfile");
my $len1=`wc -l GAUbasis`;
chomp $len1;
my @len=split(/\s+/,"$len1");
my $length=$len[0];
for (my $ig = 1; $ig <= $length; $ig +=3) {
$ele=`awk 'NR=='$ig',NR=='$ig' {print \$1}' GAUbasis`;
chomp $ele;
$here=`grep "$ele" "tfile" | wc -l`;
chomp $here;
if ($here > 0) {
my $line1=`awk 'NR=='$ig',NR=='$ig' {print \$0}' GAUbasis`;
my $line2=`awk 'NR=='$ig+1',NR=='$ig+1' {print \$0}' GAUbasis`;
my $line3=`awk 'NR=='$ig+2',NR=='$ig+2' {print \$0}' GAUbasis`;
chomp $line1;
chomp $line2;
chomp $line3;
print TMP "$line1\n";
print TMP "$line2\n";
print TMP "$line3\n";
}
}
print TMP " \n";
close TMP;
  }
 system("mv tmpfile ab.$ic1.$ic2.com");
# print IN "ab.$ic1.$ic2.com\n";
    }
close $fh;

# now the nb.* files
 system("ls -1 nb.*.0.com > nbfiles");
my $pwd=`pwd`;
chomp $pwd;
$filename='nbfiles';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 system("echo '%chk=$pwd/nb.$ic1.0.chk' > tmpfile");
# find out if there are any external charges
  $ffile[0]=0;
  if( -s "chnb.$ic1") {
  $efile=`wc -l chnb.$ic1`;
  chomp $efile;
  @ffile=split(/\s+/,$efile);
  }
# if necessary add charge to the G09 header
  system("cp IN_G09 IN_G09_temp");
  if($ffile[0] ne "0") {
   system("sed 's/nosymm/charge nosymm/' IN_G09_temp > IN_G09_current");
  }else{
   system("cp IN_G09 IN_G09_current");
  }
  system("cat IN_G09_current >> tmpfile"); 

 system("echo 'nonbonded monomer with any necessary charges\n' >> tmpfile");
 system("cat nb.$ic1.0.com >> tmpfile");
# system("echo ' ' >> tmpfile");
    if($nchgs >= 1){
  open (CHID,"chnb.$ic1");
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $m1 = ltrim($n1);
   system("cat charge.$m1.grp >> tmpfile");
   }
  close CHID;
  system("echo ' ' >> tmpfile");
    }
# now add the gen basis if necessary
  if ($ans ne 0) {
# get only the basis for the elements present
open (TMP,">>tmpfile");
my $len1=`wc -l GAUbasis`;
chomp $len1;
my @len=split(/\s+/,"$len1");
my $length=$len[0];
for (my $ig = 1; $ig <= $length; $ig +=3) {
$ele=`awk 'NR=='$ig',NR=='$ig' {print \$1}' GAUbasis`;
chomp $ele;
$here=`grep "$ele" "nb.$ic1.0.com" | wc -l`;
chomp $here;
if ($here > 0) {
my $line1=`awk 'NR=='$ig',NR=='$ig' {print \$0}' GAUbasis`;
my $line2=`awk 'NR=='$ig+1',NR=='$ig+1' {print \$0}' GAUbasis`;
my $line3=`awk 'NR=='$ig+2',NR=='$ig+2' {print \$0}' GAUbasis`;
chomp $line1;
chomp $line2;
chomp $line3;
print TMP "$line1\n";
print TMP "$line2\n";
print TMP "$line3\n";
}
}
print TMP " \n";
close TMP;

  }
 system("mv tmpfile nb.$ic1.0.com");
# print IN "nb.$ic1.0.com\n";
}
close $fh;

# now any polarisability jobs
 system("ls -1 nb.*.0-polar.com > polarfiles");
$filename='polarfiles';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 system("cp IN_POLAR tmpfile");
 system("echo 'polarisability job for nonbonded monomer\n' >> tmpfile");
 system("cat nb.$ic1.0-polar.com >> tmpfile");
# now add the gen basis if necessary
  if ($ans ne 0) {
open (TMP,">>tmpfile");
my $len1=`wc -l GAUbasis`;
chomp $len1;
my @len=split(/\s+/,"$len1");
my $length=$len[0];
for (my $ig = 1; $ig <= $length; $ig +=3) {
$ele=`awk 'NR=='$ig',NR=='$ig' {print \$1}' GAUbasis`;
chomp $ele;
$here=`grep "$ele" "nb.$ic1.0-polar.com" | wc -l`;
chomp $here;
if ($here > 0) {
my $line1=`awk 'NR=='$ig',NR=='$ig' {print \$0}' GAUbasis`;
my $line2=`awk 'NR=='$ig+1',NR=='$ig+1' {print \$0}' GAUbasis`;
my $line3=`awk 'NR=='$ig+2',NR=='$ig+2' {print \$0}' GAUbasis`;
chomp $line1;
chomp $line2;
chomp $line3;
print TMP "$line1\n";
print TMP "$line2\n";
print TMP "$line3\n";
}
}
print TMP " \n";
close TMP;
  }
 system("mv tmpfile nb.$ic1.0-polar.com");
# print IN "nb.$ic1.0-polar.com\n";
}
close $fh;
#close IN;
# end of subroutine
 }

sub SMFArunallgau {

$nfin=0;
$out="OUTLIST";
if ( -s $out ) {
$nf=`wc -l $out`;
($nfin)=split(/\s+/,$nf);
}

open(OUT,">>OUTLIST");
open(CHK,">>FCHKLIST");
open(IN,"<INLIST");
while (<IN>) {
if ( $. > $nfin ) {
$file="$_";
chomp($file);
$file=~ s/^\s+|\s+$//g;
($stem)=split(".com",$file);

system("g09 < $file > $stem.log");
$ok=`grep 'Normal termination' $stem.log | wc -l`;
 chomp($ok);
 $ok=~ s/^\s+|\s+$//g;
 if ( $ok == 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation failed for $file\n";
  close TMP;
  exit(0);
 }else{
 print OUT "$stem.log\n";
if( -s "$stem.chk") {
 system("formchk $stem.chk $stem.fchk");
 print CHK "$stem.fchk\n";
 }
 }
}
}
close OUT;
close IN;
close CHK;

$nfin=0;
$out="OUTLISTDAL";
if ( -s $out ) {
$nf=`wc -l $out`;
($nfin)=split(/\s+/,$nf);
}
open(OUT,">>OUTLISTDAL");
open(IN,"<INLISTDAL");
while (<IN>) {
if ( $. > $nfin ) {
$file="$_";
chomp($file);
($stem)=split(".mol",$file);

system("dalton -N 1 -o $stem.out $stem");
$err=`grep ERROR $stem.out | wc -l`;
chomp($err);
if ( $err > 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation failed for $stem.mol\n";
  close TMP;
  exit(0);
 }else{
 print OUT "$stem.mol\n";
 }
}
}
close OUT;
close IN;

system("rm -f OUTLIST");
system("rm -f OUTLISTDAL");
system("rm -f FCHKLIST");
}


