sub readnwcparams  {

if ( -s "xyzFILENAME" ) {
$file=`awk 'NR==1,NR==1 {print \$1}' xyzFILENAME`;
}else{
print "You must first enter the filename for the Cartesian coordinates";
system ("sleep 5");
return 1;
}
system("$EXEDIR/distinctlabels < $file ");
system("mv labels labels0");

print "                          Input variables for SMFA and NWChem\n";
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

print "\n";

print "Enter each line of NWChem input needed for the tasks,\n";
print "for example, for an SCF, DTF, MP2 or TCE-based calculation.\n";
print "Include the TASK command\n";
print "Omit the charge,\n";
print "but include 'singlet' in the SCF or DFT section.\n";
print "Begin now and end with a RETURN\n";
my $ic=0;
my $again="Y";
$nwccommands[$ic]=<STDIN>;
chomp $nwccommands[$ic];
$nwccommands[$ic]=uc($nwccommands[$ic]);
while ($again ne "N") {
my $inp=<STDIN>;
if ($inp eq "\n") {
 $again="N";
}else{
chomp $inp;
$ic=$ic+1;
$nwccommands[$ic]=uc($inp);
}
}
open (TMP1,">NWCcommands");
for (my $ie=0;$ie <= $ic;$ie++) {
print TMP1 "$nwccommands[$ie]\n";
}
close TMP1;
print "\n";
print "Does the calculation method account for electronic correlation leading to dispersion?\n";
print "(See the user's manual). Do you want to account for dispersion at long range?\n";
print "Answer Y or N\n";
$disp=<STDIN>;
chomp $disp;
$disp=uc($disp);
if($disp ne "Y") {
$disp="N";
}

#print "Enter the number of processors for each calculation:\n ";
#$nproc0=<STDIN>;
#chomp $nproc0;

print "\n";

print "NWChem allows a different basis set for each element\n";
print "Do you want all elements to have the same basis (Y or N) ? ";
$ans=<STDIN>;
chomp $ans;
$ans=uc($ans);
if($ans ne "N") {$ans="Y"};
if ($ans eq "Y") {
print "Enter the basis for all atoms\n";
$basis=<STDIN>;
chomp $basis;
open (TMP1,">NWCbasis");
print TMP1 "*   library   $basis\n";
close TMP1;
}else{
print "Enter the basis (check availability in the NWChem library)\n";
print "for the following elements\n";
open (TMP,"<labels0");
system("rm -f NWCbasis");
open (TMP1,">NWCbasis");
while (<TMP>) {
my $lab="$_";
chomp $lab;
print "$lab\n";
my $abasis=<STDIN>;
chomp $abasis;
print TMP1 "$lab","   library   ","$abasis\n";
}
close TMP;
close TMP1;
}

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
if ($ans eq "Y") {
print TMP "The basis set is\n";
print TMP "$basis\n";
}else{
print TMP "The general basis sets are\n";
close TMP;
system("cat NWCbasis >> ABSETTINGS");
      open TMP, ">>ABSETTINGS";
}
print TMP "The job type is energy (0), force (1), hessian (2), opt (3), opt=ts (4) or scan (5)\n";
       print TMP "$deriv\n";
print TMP "The NWChem commands are\n";
close TMP;
system("cat NWCcommands >> ABSETTINGS");
      open TMP, ">>ABSETTINGS";
print TMP "Is long range dispersion accounted for?\n";
       print TMP "$disp\n";
#print TMP "The number of processors for individual calculations is\n";
#       print TMP "$nproc0\n";
      close TMP;

}

sub extract_nwc {

# get the energy, force and hessian from all relevant ab initio jobs
# # run by NWChem

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

&getcoords_nwc("FRAG$ic.log");
system("cat thesecoords >> FragDerivatives");

&getenergy_nwc("FRAG$ic.log");
open(FR,">>FragDerivatives");
print FR "$TotEn\n";
close FR;
  if ($deriv >= 1) {
   &getforces_nwc("FRAG$ic.log");
   system("cat theseforces >> FragDerivatives");
  }
  if ($deriv >= 2) {
   &gethessian_nwc("FRAG$ic",$natomf);
   system("cat hessian >> FragDerivatives");
   &getdipdr_nwc("FRAG$ic.log");
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
 &getenergy_nwc("ab.$ic1.$ic2.log");
 my $TotE=$TotEn;
 &getenergy_nwc("nb.$ic1.0.log");
 $TotE=$TotE - $TotEn;
 &getenergy_nwc("nb.$ic2.0.log");
 $TotE=$TotE - $TotEn;
 my $absign=`awk '/Isg_coeff/ {print \$2}' ab.$ic1.$ic2.com`;
 chomp $absign;
 $AllTot = $AllTot + $absign * $TotE;

 if ($deriv >= 1) {
  print TMP1 "$ic1","  ","$ic2\n";
  print TMP1 "$absign\n";
  &getforces_nwc("ab.$ic1.$ic2.log");
  print TMP1 "$numberofatoms\n";
  my $nat1=$numberofatoms;
  close TMP1;
  system("cat theseforces >> abforces");
  open(TMP1,">>abforces");
  print TMP1 "$ic1\n";
  &getforces_nwc("nb.$ic1.0.log");
  print TMP1 "$numberofatoms\n";
  my $nat2=$numberofatoms;
  close TMP1;
  system("cat theseforces >> abforces");
  open(TMP1,">>abforces");
  print TMP1 "$ic2\n";
  &getforces_nwc("nb.$ic2.0.log");
    print TMP1 "$numberofatoms\n";
    my $nat3=$numberofatoms;
  close TMP1;
  system("cat theseforces >> abforces");
  open(TMP1,">>abforces");

 if ($deriv >= 2) {
  print TMP2 "$ic1","  ","$ic2\n";
  print TMP2 "$absign\n";
  &gethessian_nwc("ab.$ic1.$ic2",$nat1);
  print TMP2 "$nat1\n";
  close TMP2;
  system("cat hessian >> abhessians");
  open(TMP2,">>abhessians");
  print TMP2 "$ic1\n";
  &gethessian_nwc("nb.$ic1.0",$nat2);
  print TMP2 "$nat2\n";
  close TMP2;
  system("cat hessian >> abhessians");
  open(TMP2,">>abhessians");
  print TMP2 "$ic2\n";
  &gethessian_nwc("nb.$ic2.0",$nat3);
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

sub getcoords_nwc {
$filename="$_[0]";
 open(FRAGJOB,"<$filename");
 open(TMP,">thesecoords");
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
open(TMP,">thesedipdr");
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
 system("grep 'Total SCF energy' $filename > junk");
 }
 $dft=`grep 'Total DFT energy' $filename | wc -l`;
 chomp($dft);
 $dft=~ s/^\s+|\s+$//g;
 if( $dft > 0 ) {
 system("grep 'Total DFT energy' $filename > junk");
 }
 $line1=`awk 'NR==1,NR==1 {print \$0}' junk`;
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

}

sub getforces_nwc {
$numberofatoms=0;
 $skip=0;
 $writit=0;
 open(TMP,">theseforces");
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
system("cp $file1.hess hessian");
#@file=split(/\./,"$file1");
#system("cp $file[0].hess hessian");
}


sub makenbgdma_nwc {
# actually we cannot use gdma for NWChem
# so we use the ESP charges only

my $nfrags=`ls nb.*.0.nw | wc -l`;
 for (my $ic=1; $ic<=$nfrags; $ic++)
  {
 open(FH1,"<nb.$ic.0.nw");
 my $skip=0;
 my $writit=0;
 $icount=0;
 while (<FH1>) {
   if (/symmetry c1/) {
    $skip=$. + 1;
   }
 if ($. == $skip) {$writit=1};
 if(/END/)  {$writit=0};
      if( $writit == 1) {
         $line1="$_";
         chomp $line1;
         @line=split(/\s+/,$line1);
         $icount=$icount+1;
         $x[$icount]=$line[1];     
         $y[$icount]=$line[2];     
         $z[$icount]=$line[3];     
      }
  } 
close FH1;
 $icount=0;
 $natom=`awk 'NR==1,NR==1 {print \$1}' nb.$ic.0.esp`;
 chomp $natom;
 system("tail -$natom nb.$ic.0.esp > tfile");
 open(FH1,"<tfile");
 while (<FH1>) {
  $line1="$_";
  chomp $line1;
  @line=split(/\s+/,$line1);
  $icount=$icount+1;
  $lab[$icount]=$line[0];
  $q[$icount]=$line[4];
 }
close FH1;

open(TMP,">nb.$ic.0.cart");
print TMP "  The distributed Cartesian multipoles\n";
print TMP "$natom\n";
for (my $ie=1;$ie <= $natom;$ie++) {
 print TMP "$lab[$ie]\n";
 print TMP "$x[$ie]\n";
 print TMP "$y[$ie]\n";
 print TMP "$z[$ie]\n";
 print TMP "Charge\n";
 print TMP "$q[$ie]\n";
 print TMP "Dipole\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "Quodrupole\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "Octapole\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
print TMP "Hexadecapole\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
 print TMP "0\n";
}
close TMP;
}

}



sub SMFArunnwc {

  my $nfrgs = `ls COORD* | wc -l`;
  chomp $nfrgs;
 for (my $ic=1; $ic<=$nfrgs; $ic++)
  {
  system("nwchem FRAG$ic.nw > FRAG$ic.log");
 $ok=`grep 'CITATION' FRAG$ic.log | wc -l`;
 chomp($ok);
 $ok=~ s/^\s+|\s+$//g;
 if ( $ok == 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation failed for FRAG$ic\n";
  close TMP;
  exit(0);
 }
  }

 system("ls -1 ab*.nw > abfiles");
$filename='abfiles';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 my $ic2=$fields[2];
 system("nwchem ab.$ic1.$ic2.nw > ab.$ic1.$ic2.log");
 $ok=`grep 'CITATION' ab.$ic1.$ic2.log | wc -l`;
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

 system("ls -1 nb*.0.nw > nbfiles");
$filename='nbfiles';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 system("nwchem nb.$ic1.0.nw > nb.$ic1.0.log");
 $ok=`grep 'CITATION' nb.$ic1.0.log | wc -l`;
 chomp($ok);
 $ok=~ s/^\s+|\s+$//g;
 if ( $ok == 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation failed for nb.$ic1.0\n";
  close TMP;
  exit(0);
 }
 }
close $fh;

# the polarization jobs are run using DALTON
# by the script rundalpolar
}

sub SMFAnwcinputs_MAC {

    $nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;

$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
my $n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
$n++;
}
$n=$n-1;
#$nproc0=$store[$n];
close $fh;

for (my $ic=0; $ic<=$n; $ic++) {
 if ("$store[$ic]" eq "The electronic structure method is") {
    my $ie=$ic+1;
    $method=$store[$ie];
 }
}

#for (my $ic=0; $ic<=$n; $ic++) {
# if ("$store[$ic]" eq "The job type is energy \(0\), force \(1\) or hessian \(2\)") {
#    $ie=$ic+1;
#    $deriv=$store[$ie];
# }
#}

$deriv =`awk 'NR==2,NR==2 {print \$1}' IN_JOBTYPE`;
chomp $deriv;

#if ($deriv == 0) {my $prop = "energy"};
#if ($deriv == 1) {my $prop = "gradient"};
#if ($deriv == 2) {my $prop = "hessian"};

#open(IN,">INLIST");

# start with the main fragments
  my $nfrgs = `ls COORD* | wc -l`;
  chomp $nfrgs;
 for (my $ic=1; $ic<=$nfrgs; $ic++)
  {
  open(TMP,">tmpfile");
  my $my_coord = "COORD${ic}";
  print TMP "START FRAG$ic\n";
  print TMP "ECHO\n";
  $ch=`awk 'NR==1,NR==1{print \$1}' COORD$ic`;
  $mult=`awk 'NR==1,NR==1{print \$3}' COORD$ic`;
  chomp $ch;
  chomp $mult;
  if ("$mult" eq "1") {
   $mform="SINGLET";
  }
  if ("$mult" eq "2") {
   $mform="doublet \;  uhf";
  }
  print TMP "charge $ch\n";
  print TMP "geometry nocenter\n";
  print TMP "symmetry c1\n";
  my $lf=`wc -l COORD$ic`;
  chomp $lf;
  my @kf=split(/\s/,"$lf");
  $mf=$kf[0];
  $mf=$mf-1;
  close TMP;
  system("tail -$mf COORD$ic >> tmpfile");
  open(TMP,">>tmpfile");
  print TMP "END\n";
if($nchgs >= 1){
  open (CHID,"chFRAG.$ic");
  print TMP "bq\n";
  close TMP; 
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $m1 = &ltrim($n1);
   system("cat charge.$m1.grp >> tmpfile");
   }
  close CHID;
  open(TMP,">>tmpfile");
    }
  if($nchgs >= 1){print TMP "END\n"};

  print TMP "BASIS\n";
  close TMP;
  system("cat NWCbasis >> tmpfile");
  open(TMP,">>tmpfile");
  print TMP "END\n";
 open(FH,"<NWCcommands");
 while (<FH>) {
 chomp($_);
  if( $mult == 2 ) {
   s/SINGLET/$mform/;
  }  
  print TMP "$_\n";
 }
  close TMP;
  system("mv tmpfile FRAG$ic.nw");
#  print IN "FRAG$ic.nw\n";
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
  open(TMP,">tmpfile");
  print TMP "START ab.$ic1.$ic2\n";
  print TMP "ECHO\n";
# get the sign of the ab nb combination
  $_=`head -1 ab.$ic1.$ic2.com`;
  s/!Isg_coeff=/\$comment!Isg_coeff=/;
  $comm=$_;
  chomp $comm;
  print TMP "#$comm\n";
  $ch=`awk 'NR==3,NR==3{print \$1}' ab.$ic1.$ic2.com`;
  $mult=`awk 'NR==3,NR==3{print \$2}' ab.$ic1.$ic2.com`;
  chomp $ch;
  chomp $mult;
  if ($mult == 1) {$mform="SINGLET"};
  if ($mult == 2) {$mform="doublet ;  uhf"};
  print TMP "charge $ch\n";
  print TMP "geometry nocenter\n";
  print TMP "symmetry c1\n";
  close TMP;
 my $lf=`wc -l ab.$ic1.$ic2.com`;
 chomp $lf;
 my @kf=split(/\s/,"$lf");
 my $mf=$kf[0];
 $mf=$mf-1;
 system("head -$mf ab.$ic1.$ic2.com > tfile");
 $mf=$mf-3;
 system("tail -$mf tfile  >> tmpfile");
 open(TMP,">>tmpfile");
 print TMP "END\n";
if($nchgs >= 1){
  open (CHID,"chab.$ic1.$ic2");
  print TMP "bq\n";
  close TMP;
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $m1 = ltrim($n1);
   system("cat charge.$m1.grp >> tmpfile");
   }
  close CHID;
  open(TMP,">>tmpfile");
    }
  if($nchgs >= 1){print TMP "END\n"};
  print TMP "BASIS\n";
  close TMP;
  system("cat NWCbasis >> tmpfile");
  open(TMP,">>tmpfile");
  print TMP "END\n";
 open(FH,"<NWCcommands");
 while (<FH>) {
 chomp($_);
  if( $mult == 2 ) {
   s/SINGLET/$mform/;
  }
  print TMP "$_\n";
 }
  close TMP;
  system("mv tmpfile ab.$ic1.$ic2.nw");
#  print IN "ab.$ic1.$ic2.nw\n";
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
  open(TMP,">tmpfile");
  print TMP "START nb.$ic1.0\n";
  print TMP "ECHO\n";
  $ch=`awk 'NR==1,NR==1{print \$1}' nb.$ic1.0.com`;
  $mult=`awk 'NR==1,NR==1{print \$2}' nb.$ic1.0.com`;
  chomp $ch;
  chomp $mult;
  if ($mult == 1) {$mform="SINGLET"};
  if ($mult == 2) {$mform="doublet ;  uhf"};
  print TMP "charge $ch\n";
  print TMP "geometry nocenter\n";
  print TMP "symmetry c1\n";
  close TMP;
 my $lf=`wc -l nb.$ic1.0.com`;
 chomp $lf;
 my @kf=split(/\s/,"$lf");
 my $mf=$kf[0];
 $mf=$mf-1;
 system("head -$mf nb.$ic1.0.com > tfile");
 $mf=$mf-1;
 system("tail -$mf tfile >> tmpfile");
 open(TMP,">>tmpfile");
 print TMP "END\n";
if($nchgs >= 1){
  open (CHID,"chnb.$ic1");
  print TMP "bq\n";
  close TMP;
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $m1 = ltrim($n1);
   system("cat charge.$m1.grp >> tmpfile");
   }
  close CHID;
  open(TMP,">>tmpfile");
    }
  if($nchgs >= 1){print TMP "END\n"};
  print TMP "BASIS\n";
  close TMP;
  system("cat NWCbasis >> tmpfile");
  open(TMP,">>tmpfile");
  print TMP "END\n";
 open(FH,"<NWCcommands");
 while (<FH>) {
 chomp($_);
  if( $mult == 2 ) {
   s/SINGLET/$mform/;
  }
  print TMP "$_\n";
 }
  print TMP "ESP\n";
  print TMP "recalculate\n";
  print TMP "END\n";
  print TMP "TASK ESP";
  close TMP;
  system("mv tmpfile nb.$ic1.0.nw");
#  print IN "nb.$ic1.0.nw\n";
}
close $fh;

# the nb polar files are constructed for DALTON
# as NWChem does not have polarizabilities
# see rundalpolar
#close IN;
}

sub Lev0_chg_MAC_iter_nwc {

system("clear");

$count=1;
$iter=3;

while ($count <= $iter) {

&SMFAnwcnpa_MAC;
&extractch_nwc;
      $count++;
}
      system("rm charge*.db");
      system("rm charge*.plt");
      return 1;
}

sub extractch_nwc  {

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

system("tail -$Nat charge.$ic.esp > tmpfile");
      open (NPA_IN,"<tmpfile");
      open (NPA_OUT,">charge.$ic.npa");
      print NPA_OUT "charge.$ic.coord\n";
      print NPA_OUT "$Nat\n";
      $n=0;
      while(<NPA_IN>) {
      $line="$_";
      @bits=split(/\s+/,"$line");
      $q=$bits[4];
      print NPA_OUT "$x[$n] "," $y[$n] "," $z[$n] "," $q\n";
      $n++;
}
      close NPA_OUT;
      system("$EXEDIR/adjustcaponly < charge.$ic.npa");
}
# end sub
}

sub SMFAnwcnpa_MAC {

# Although the name of this subroutine contains npa
# the npa charges are not available on nwc, so we use
# the ESP module of nwc

&SMFAmkchinputs_nwc;
&runchseq_nwc;

}

sub SMFAmkchinputs_nwc {
 my @jobs = `ls charge.*.coord`;
 foreach my $chg_coord (@jobs) {
  @fields=split(/\./,$chg_coord);
  $ic=$fields[1];

  open(TMP,">tmpfile");
  my $my_coord = "charge.$ic.coord";
  print TMP "START charge.$ic\n";
  print TMP "ECHO\n";
  $ch=`awk 'NR==1,NR==1{print \$1}' charge.$ic.coord`;
  $mult=`awk 'NR==1,NR==1{print \$2}' charge.$ic.coord`;
  chomp $ch;
  chomp $mult;
  if ("$mult" eq "1") {
   $mform="SINGLET";
  }
  if ("$mult" eq "2") {
   $mform="doublet \;  uhf";
  }
  print TMP "charge $ch\n";
  print TMP "geometry nocenter\n";
  print TMP "symmetry c1\n";
  my $lf=`wc -l charge.$ic.coord`;
  chomp $lf;
  my @kf=split(/\s/,"$lf");
  $mf=$kf[0];
  $mf=$mf-1;
  close TMP;
  system("tail -$mf charge.$ic.coord >> tmpfile");
  open(TMP,">>tmpfile");
  print TMP "END\n";
  print TMP "BASIS\n";
  close TMP;
  system("cat NWCbasis >> tmpfile");
  open(TMP,">>tmpfile");
  print TMP "END\n";
  print TMP "SCF\n";
  print TMP " direct\n";
  print TMP " $mform\n";
  print TMP "END\n";
if($nchgs > 1){
  open (CHID,"chL0.$ic");
  print TMP "bq\n";
  close TMP;
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $m1 = ltrim($n1);
   system("cat charge.$m1.grp >> tmpfile");
   }
  close CHID;
  open(TMP,">>tmpfile");
  print TMP "END\n";
    }
  print TMP "ESP\n";
  print TMP "recalculate\n";
  print TMP "END\n";
  print TMP "TASK SCF ENERGY\n";
  print TMP "TASK ESP";
  close TMP;
  system("mv tmpfile charge.$ic.nw");
}

# end sub
}

sub runchseq_nwc  {
    $nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
    chomp $nchgs;

 my @jobs = `ls charge.*.coord`;
 foreach my $chg_coord (@jobs) {
  @fields=split(/\./,$chg_coord);
  $ic=$fields[1];
  system("mpirun nwchem charge.$ic.nw > charge.$ic.out");
 $ok=`grep 'CITATION' charge.$ic.out | wc -l`;
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


sub SMFArunallnwc {

$nfin=0;
$out="OUTLIST";
if ( -s $out ) {
$nf=`wc -l $out`;
($nfin)=split(/\s+/,$nf);
}

open(OUT,">>OUTLIST");
open(IN,"<INLIST");
while (<IN>) {
if ( $. > $nfin ) {
$file="$_";
chomp($file);
$file=~ s/^\s+|\s+$//g;
($stem)=split(".com",$file);
$file=join("",$stem,".nw");

system("mpirun nwchem $file > $stem.log");
$ok=`grep 'CITATION' $stem.log | wc -l`;
 chomp($ok);
 $ok=~ s/^\s+|\s+$//g;
 if ( $ok == 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation failed for $file\n";
  close TMP;
  exit(0);
 }else{
 print OUT "$stem.log\n";
 }
}
}
close OUT;
close IN;

#$nfin=0;
#$out="OUTLISTDAL";
#if ( -s $out ) {
#$nf=`wc -l $out`;
#($nfin)=split(/\s+/,$nf);
#}
#open(OUT,">>OUTLISTDAL");
#open(IN,"<INLISTDAL");
#while (<IN>) {
#if ( $. > $nfin ) {
#$file="$_";
#chomp($file);
#($stem)=split(".mol",$file);
#system ("ln -s IN_DISP $stem.dal");
#system("dalton -N 1 -o $stem.out $stem");
#$err=`grep ERROR $stem.out | wc -l`;
#chomp($err);
#if ( $err > 0 ) {
#  open(TMP,">>OUT_SMFA");
#  print TMP "calculation failed for $stem.mol\n";
#  close TMP;
#  exit(0);
# }else{
# print OUT "$stem.mol\n";
# }
#}
#}
#close OUT;
#close IN;

system("rm -f OUTLIST");
#system("rm -f OUTLISTDAL");
}















