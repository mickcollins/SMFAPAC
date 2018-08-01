sub readgamparams {

if ( -s "xyzFILENAME" ) {
$file=`awk 'NR==1,NR==1 {print \$1}' xyzFILENAME`;
}else{
print "You must first enter the filename for the Cartesian coordinates\n";
system ("sleep 5");
return 1;
}
system("$EXEDIR/distinctlabels < $file ");
system("mv labels labels0");

print "                              GAMESS input\n";
print "\n";
$deriv=0;
print "Enter the calculation type as an integer:\n";
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


      open TMP, ">IN_JOBTYPE";
      print TMP "Enter 0, 1, 2, 3, 4 for energy, gradient, hessian, opt or TS\n";
      print TMP "$deriv\n";
      close TMP;


print "\n";
print "               \$CONTRL group variables\n";
print "\n";
print "Here you enter necessary parts of the \$CONTRL group, ONE item per line\n";
print "For example\n";
print "SCFTYP=RHF\n";
print "MPLEVL=2\n";
print "\n";
print "Omit the RUNTYP variable\n";
print "Omit the charge \(ICHARG\) and multiplicity \(MULT\)\n";
print "Omit \$BASIS\n";
print "Omit COORD as COORD=UNIQUE is used by SMFA\n";
print "Omit the \$DATA section\n";
print "The format MUST BE EXACTLY as GAMESS requires\n";
print "Begin now and end with a RETURN\n";
my $ic=0;
my $again="Y";
$gamcommands[$ic]=<STDIN>;
chomp $gamcommands[$ic];
$gamcommands[$ic]=uc($gamcommands[$ic]);
while ($again ne "N") {
my $inp=<STDIN>;
if ($inp eq "\n") {
 $again="N";
}else{
chomp $inp;
$ic=$ic+1;
$gamcommands[$ic]=uc($inp);
}
}
open (TMP1,">GAMcommands");
for (my $ie=0;$ie <= $ic;$ie++) {
print TMP1 "$gamcommands[$ie]\n";
}
close TMP1;
print "\n";
print "                 Additional GAMESS variables\n";
print "\n";
print "Here you enter other GAMESS variable groups (if any), eg \$SYSTEM, \$SCF, etc\n";
print "The format MUST BE EXACTLY as GAMESS requires (don't leave out compulsory spaces)\n";
print "Begin now and end with a RETURN\n";
my $ic=0;
my $again="Y";
$gamcommands[$ic]=<STDIN>;
chomp $gamcommands[$ic];
$gamcommands[$ic]=uc($gamcommands[$ic]);
while ($again ne "N") {
my $inp=<STDIN>;
if ($inp eq "\n") {
 $again="N";
}else{
chomp $inp;
$ic=$ic+1;
$gamcommands[$ic]=uc($inp);
}
}
open (TMP1,">GAMextras");
for (my $ie=0;$ie <= $ic;$ie++) {
print TMP1 "$gamcommands[$ie]\n";
}
close TMP1;

print "\n";
print "Does the chosen electronic structure  method specified above account for\n";
print "electronic correlation leading to dispersion (see the user's manual)?\n";
print "Do you therefore want to account for dispersion at long range?\n";
print "Answer Y or N\n";
$disp=<STDIN>;
chomp $disp;
$disp=uc($disp);
if($disp ne "Y") {
$disp="N";
}

print "\n";
system("rm -f GAMbasis");
print "GAMESS allows a different basis set for each chemical element\n";
print "Do you want all elements to have the same basis (Y or N) ? ";
my $ans=<STDIN>;
chomp $ans;
$ans=uc($ans);
if($ans ne "N") {$ans="Y"};
print "\n";
open (TMP1,">GAMbasis");
print TMP1 "$ans\n";
if ($ans eq "Y") {
print "Enter the common basis for all atoms\n";
print "You MUST follow the format required by GAMESS\n";
print "in 'card sequence U' of the \$DATA group\n";
print "Enter one or more lines and finish with a (blank) RETURN\n";
my $again="Y";
while ($again eq "Y") {
$basis=<STDIN>;
if ($basis eq "\n") {
 $again="N";
}else{
chomp $basis;
print TMP1 "$basis\n";
}
}
print TMP1 "\n";
}else{
print "To enter the basis for each chemical element,\n";
print "you MUST follow the format required by GAMESS";
print "in 'card sequence U' of the \$DATA group\n";
print "and finish each basis set with a (blank) RETURN\n";
print "Enter the basis (check availability in the GAMESS library)\n";
print "for the following elements\n";
open (TMP,"<labels0");
while (<TMP>) {
my $lab="$_";
chomp $lab;
print "$lab\n";
print TMP1 "$lab\n";
my $again="Y";
while ($again eq "Y") {
$basis=<STDIN>;
if ($basis eq "\n") {
 $again="N";
}else{
chomp $basis;
print TMP1 "$basis\n";
}
}
print TMP1 "\n";
}
close TMP;
close TMP1;
}
# get the basis set used for the polarizability in DALTON
print "\n";
print "SMFA often needs to calculate the static dipole polarizabilty of each functional group.\n";
print "SMFA uses the DALTON program, rather than GAMESS to do this.\n";
print "The 6-311+G(d,p) basis set is considered adequate for this purpose, and is the default.\n";
print "If you want to over-ride this default, you must enter a basis set in a format that\n";
print "is recognsed by DALTON.\n";
print "Enter the chosen basis set, or hit the RETURN key to accept the default:\n";
$dalbasis="6-311++G**";
$basisin=<STDIN>;
if ($basisin eq "\n") {
}else{
chomp $basisin;
$basisin=~ s/^\s+|\s+$//g;
$dalbasis=$basisin;
}


if($deriv==3) {&optinput};
if($deriv==4) {&optinput};
if($deriv==4) {&TSinput};
if($deriv==5) {&scaninput};

#print "\n";
#print "Enter the number of processors for each calculation:\n ";
#$nproc0=<STDIN>;
#chomp $nproc0;

# Note that the $contrl card with charge and multiplicity is missing

      open TMP, ">ABSETTINGS";
       print TMP "The quantum chemistry program package is\n";
       print TMP "$package\n";
       print TMP "The job type is energy \(0\), force \(1\), hessian \(2\), opt \(3\), opt=TS \(4\), scan \(5\)\n";
       print TMP "$deriv\n";
       print TMP "The GAMESS  'cards' are\n";
      close TMP;
       system("cat GAMcommands >> ABSETTINGS");
       system("cat GAMextras >> ABSETTINGS");
      open TMP, ">>ABSETTINGS"; 
       print TMP "The basis set(s)\n";
      close TMP;
my $len1=`wc -l GAMbasis`;
my @len=split(/\s+/,"$len1");
my $length=$len[0];
$length=$length-1;
       system("tail -$length GAMbasis >> ABSETTINGS");
      open TMP, ">>ABSETTINGS";
       print TMP "Is long range dispersion accounted for?\n";
       print TMP "$disp\n";
       print TMP "The Dalton basis for polarizabilities is\n";
       print TMP "$dalbasis\n";
#       print TMP "The number of processors for individual calculations is\n";
#       print TMP "$nproc0\n";
       close TMP;

}


sub getStonegam {
$stem="@_";
$filename=join("",$stem,"log");
$outfile=join("",$stem,"Stone");
$cartfile=join("",$stem,"cart");

open(TMP,"<$filename");
open(STO,">$outfile");

$skip=0;
$writit=0;
$ic=0;
while (<TMP>) {
if(/NET CHARGES AT POINTS                    COORDINATES/) {
$skip= $.+ 3;
}
if($skip == $.) {$writit=1};
if(/THE DISTRIBUTED MULTIPOLE ANALYSIS IS/) {$writit=0}
if($writit == 1) {
if($_ eq "\n") {
}else{
 $ic=$ic+1;
 chomp $_;
 $lines[$ic]=$_;
}
}
}
close TMP;
print STO "$ic\n";
for ($ie=1;$ie <= $ic;$ie++) {print STO "$lines[$ie]\n"};

open(TMP,"<$filename");
$skip=0;
$writit=0;
$ic=0;
while (<TMP>) {
if(/FIRST MOMENTS AT POINTS/) {
$skip= $.+ 3;
}
if($skip == $.) {$writit=1};
if(/SECOND MOMENTS AT POINTS/) {$writit=0}
if($writit == 1) {
if($_ eq "\n") {
}else{
 $ic=$ic+1;
 chomp $_;
 $lines[$ic]=$_;
}
}
}
close TMP;
for ($ie=1;$ie <= $ic;$ie++) {print STO "$lines[$ie]\n"};

open(TMP,"<$filename");
$skip=0;
$writit=0;
$ic=0;
while (<TMP>) {
if(/SECOND MOMENTS AT POINTS/) {
$skip= $.+ 3;
}
if($skip == $.) {$writit=1};
if(/THIRD MOMENTS AT POINTS/) {$writit=0}
if($writit == 1) {
if($_ eq "\n") {
}else{
 $ic=$ic+1;
 chomp $_;
 $lines[$ic]=$_;
}
}
}
close TMP;
for ($ie=1;$ie <= $ic;$ie++) {print STO "$lines[$ie]\n"};

open(TMP,"<$filename");
$skip=0;
$writit=0;
$ic=0;
while (<TMP>) {
if(/THIRD MOMENTS AT POINTS/) {
$skip= $.+ 4;
}
if($skip == $.) {$writit=1};
if(/STEP CPU TIME/) {$writit=0}
if($writit == 1) {
if($_ eq "\n") {
}else{
 $ic=$ic+1;
 chomp $_;
 $lines[$ic]=$_;
}
}
}
close TMP;
for ($ie=1;$ie <= $ic;$ie++) {print STO "$lines[$ie]\n"};

close STO;

system("$EXEDIR/gamstone2cart < $outfile > $cartfile");

}

sub extract_gam {

# get the energy, force and hessian from all re;evant ab initio jobs
# run by GAMESS

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

&getcoords_gam("FRAG$ic.log");
system("cat thesecoords >> FragDerivatives");

&getenergy_gam("FRAG$ic.log");
open(FR,">>FragDerivatives");
print FR "$TotEn\n";
close FR;
  if ($deriv >= 1) {
   &getforces_gam("FRAG$ic.dat");
   system("cat theseforces >> FragDerivatives");
  }
  if ($deriv >= 2) {
   &gethessian_gam("FRAG$ic.log","$numberofatoms");
   system("cat hessian >> FragDerivatives");
   &getdipdr_gam("FRAG$ic.log");
   system("echo '$numberofatoms' >> FragDipderivs");
   system("cat thesedipdr >> FragDipderivs");
  }
 }
# now the ab initio nonbonded jobs
system("ls -1 ab.*.dat > abfiles");
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
 &getenergy_gam("ab.$ic1.$ic2.log");
 my $TotE=$TotEn;
 &getenergy_gam("nb.$ic1.0.log");
 $TotE=$TotE - $TotEn;
 &getenergy_gam("nb.$ic2.0.log");
 $TotE=$TotE - $TotEn;
 my $absign=`awk '/Isg_coeff/ {print \$2}' ab.$ic1.$ic2.com`;
 chomp $absign;
 $AllTot = $AllTot + $absign * $TotE;

 if ($deriv >= 1) {
  print TMP1 "$ic1","  ","$ic2\n";
  print TMP1 "$absign\n";
  &getforces_gam("ab.$ic1.$ic2.dat");
  print TMP1 "$numberofatoms\n";
  my $nat1=$numberofatoms;
  close TMP1;
  system("cat theseforces >> abforces");
  open(TMP1,">>abforces");
  print TMP1 "$ic1\n";
  &getforces_gam("nb.$ic1.0.dat");
  print TMP1 "$numberofatoms\n";
  my $nat2=$numberofatoms;
  close TMP1;
  system("cat theseforces >> abforces");
  open(TMP1,">>abforces");
  print TMP1 "$ic2\n";
  &getforces_gam("nb.$ic2.0.dat");
    print TMP1 "$numberofatoms\n";
    my $nat3=$numberofatoms;
  close TMP1;
  system("cat theseforces >> abforces");
  open(TMP1,">>abforces");

 if ($deriv >= 2) {
  print TMP2 "$ic1","  ","$ic2\n";
  print TMP2 "$absign\n";
  &gethessian_gam("ab.$ic1.$ic2.log",$nat1);
  print TMP2 "$nat1\n";
  close TMP2;
  system("cat hessian >> abhessians");
  open(TMP2,">>abhessians");
  print TMP2 "$ic1\n";
  &gethessian_gam("nb.$ic1.0.log",$nat2);
  print TMP2 "$nat2\n";
  close TMP2;
  system("cat hessian >> abhessians");
  open(TMP2,">>abhessians");
  print TMP2 "$ic2\n";
  &gethessian_gam("nb.$ic2.0.log",$nat3);
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

sub getcoords_gam {

 open(FRAGJOB,"<$_[0]");
 open(TMP,">thesecoords");

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
open(TMP,">thesedipdr");
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
 open(TMP,">theseforces");
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
 open(TMP,">tfile");
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

system("$EXEDIR/convertgamhessian < tfile > hessian");

}

sub makenbgdma_gam {

my $nfrgs = `ls nb.*.0.log | wc -l`;
chomp $nfrgs;

for (my $m=1;$m <= $nfrgs;$m++) {
$stem="nb.$m.0.";
$filename=join("",$stem,"log");
$outfile=join("",$stem,"Stone");
$cartfile=join("",$stem,"cart");



open(TMP,"<$filename");
open(STO,">$outfile");

$skip=0;
$writit=0;
$ic=0;
$once=0;
while (<TMP>) {
if(/NET CHARGES AT POINTS                    COORDINATES/) {
$skip= $.+ 3;
}
if($skip == $.) {$writit=1};
if(/THE DISTRIBUTED MULTIPOLE ANALYSIS IS/) {
 $writit=0;
 $once=1;
}
if($writit == 1 && $once == 0) {
if($_ eq "\n") {
}else{
 $ic=$ic+1;
 chomp $_;
 $lines[$ic]=$_;
}
}
}
close TMP;
print STO "$ic\n";
for ($ie=1;$ie <= $ic;$ie++) {print STO "$lines[$ie]\n"};

open(TMP,"<$filename");
$skip=0;
$writit=0;
$ic=0;
$once=0;
while (<TMP>) {
if(/FIRST MOMENTS AT POINTS/) {
$skip= $.+ 3;
}
if($skip == $.) {$writit=1};
if(/SECOND MOMENTS AT POINTS/) {
 $writit=0;
 $once=1;
}
if($writit == 1 && $once == 0) {
if($_ eq "\n") {
}else{
 $ic=$ic+1;
 chomp $_;
 $lines[$ic]=$_;
}
}
}
close TMP;
for ($ie=1;$ie <= $ic;$ie++) {print STO "$lines[$ie]\n"};

open(TMP,"<$filename");
$skip=0;
$writit=0;
$ic=0;
$once=0;
while (<TMP>) {
if(/SECOND MOMENTS AT POINTS/) {
$skip= $.+ 3;
}
if($skip == $.) {$writit=1};
if(/THIRD MOMENTS AT POINTS/) {
 $writit=0;
 $once=1;
}
if($writit == 1 && $once == 0) {
if($_ eq "\n") {
}else{
 $ic=$ic+1;
 chomp $_;
 $lines[$ic]=$_;
}
}
}
close TMP;
for ($ie=1;$ie <= $ic;$ie++) {print STO "$lines[$ie]\n"};

open(TMP,"<$filename");
$skip=0;
$writit=0;
$ic=0;
$once=0;
while (<TMP>) {
if(/THIRD MOMENTS AT POINTS/) {
$skip= $.+ 4;
}
if($skip == $.) {$writit=1};
if(/STEP CPU TIME/) {
 if($writit == 1) {
  $writit=0;
  $once=1;
 }
}
if($writit == 1 && $once == 0) {
if($_ eq "\n") {
}else{
 $ic=$ic+1;
 chomp $_;
 $lines[$ic]=$_;
}
}
}
close TMP;
for ($ie=1;$ie <= $ic;$ie++) {print STO "$lines[$ie]\n"};

close STO;

system("$EXEDIR/gamstone2cart < $outfile > $cartfile");

}
}



sub SMFArungam {

  my $nfrgs = `ls COORD* | wc -l`;
  chomp $nfrgs;
 for (my $ic=1; $ic<=$nfrgs; $ic++)
  {
  system("rungms FRAG$ic.inp \$PBS_NCPUS > FRAG$ic.log");
 $ok=`grep 'TERMINATED NORMALLY' FRAG$ic.log | wc -l`;
 chomp($ok);
 $ok=~ s/^\s+|\s+$//g;
 if ( $ok == 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation failed for FRAG$ic\n";
  close TMP;
  exit(0);
 }

  }

 system("ls -1 ab*.inp > abfiles");
$filename='abfiles';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 my $ic2=$fields[2];
 system("rungms ab.$ic1.$ic2.inp \$PBS_NCPUS > ab.$ic1.$ic2.log");
 $ok=`grep 'TERMINATED NORMALLY' ab.$ic1.$ic2.log | wc -l`;
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

 system("ls -1 nb*.0.inp > nbfiles");
$filename='nbfiles';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 system("rungms nb.$ic1.0.inp \$PBS_NCPUS > nb.$ic1.0.log");
 $ok=`grep 'TERMINATED NORMALLY' nb.$ic1.0.log | wc -l`;
 chomp($ok);
 $ok=~ s/^\s+|\s+$//g;
 if ( $ok == 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation failed for nb.$ic1.0\n";
  close TMP;
  exit(0);
 }
# system("formchk nb.$ic1.0.chk nb.$ic1.0.fchk");
 }
close $fh;

# the polarization jobs are run using DALTON
# by the script rundalpolar
}

sub SMFAgaminputs_MAC {

    $nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;
# nchgs = 0 for GAMESS
$nchgs = 0;

$deriv=`awk 'NR==2,NR==2 {print \$1}' IN_JOBTYPE`;
chomp $deriv;
if ($deriv == 0) {$prop = "RUNTYP=ENERGY"};
if ($deriv == 1) {$prop = "RUNTYP=GRADIENT"};
if ($deriv == 2) {$prop = "RUNTYP=HESSIAN"};

# get the original value of deriv, as this may be opt,ts or scan
$oldderiv=`awk 'NR==2,NR==2 {print \$1}' IN_JOBTYPE_original`;
chomp $oldderiv;

# get CONTRL group commands
open(TMP,"<GAMcommands");
$ic=0;
while(<TMP>) {
my $line="$_";
chomp $line;
$ic=$ic+1;
$prop=join(" ",$prop,$line);
}
open(JUNK,">gamjunk");
print JUNK "$prop\n";
close JUNK;
system("sed 's/RHF/UHF/' gamjunk > gamjunk2");
open(JUNK2,"<gamjunk2");
while (<JUNK2>) {
$prop2="$_";
}
close JUNK2;
chomp($prop2);

  $ans=`awk 'NR==1,NR==1 {print \$1}' GAMbasis`;
  chomp $ans;
  my $lf=`wc -l GAMbasis`;
  chomp $lf;
  my @kf=split(/\s/,"$lf");
  $mf=$kf[0];
  $mf=$mf-1;
  system("tail -$mf GAMbasis > gbfile");

#open(IN,">INLIST");

# start with the main fragments
  my $nfrgs = `ls COORD* | wc -l`;
  chomp $nfrgs;
 for (my $ic=1; $ic<=$nfrgs; $ic++)
  {
  $ch=`awk 'NR==1,NR==1{print \$1}' COORD$ic`;
  $mult=`awk 'NR==1,NR==1{print \$3}' COORD$ic`;
  chomp $ch;
  chomp $mult;
  $mult=~ s/^\s+|\s+$//g;
  open(TMP,">tmpfile");
  if($mult == 1) {
  print TMP " \$CONTRL $prop MULT=$mult ICHARG=$ch \$END\n";
  }
  if($mult == 2) {
  print TMP " \$CONTRL $prop2 MULT=$mult ICHARG=$ch \$END\n";
  }
  close TMP;

  if ( -s "GAMextras") {
   if( $deriv == 2 && $oldderiv > 2 ) {
# this is the hessian at the start of an opt or ts
# so we only use scf for gamess
    system("grep SYSTEM GAMextras >> tmpfile");
   }else{
   system("cat GAMextras >> tmpfile");
   }
  }
  open(TMP,">>tmpfile");
  print TMP " \$DATA\n";
  print TMP " GAMESS job\n";
  print TMP "c1\n";
  
  my $lf=`wc -l COORD$ic`;
  chomp $lf;
  my @kf=split(/\s/,"$lf");
  $mf=$kf[0];
  $mf=$mf-1;
  system("tail -$mf COORD$ic > tfile");
  system("$EXEDIR/convertcoords4games < tfile > gamesscoords");
# have to interleave gamesscoords with atombasis and blanl line
  open(TMP1,"<gamesscoords");

  while (<TMP1>) {
   my $line1=$_;
   chomp $line1;
   my @bits=split(/\s+/,$line1);
   $ele=$bits[0];
   print TMP "$line1\n";
if ($ans eq "N") {
  open(GAMB,"<gbfile");
  $writit=0;
  $skip=0;
  while (<GAMB>) {
   my $line2=$_;
   chomp $line2;
   if($line2 eq "") {
    $writit=0;
   }else{
   my @bits=split(/\s+/,$line2);
   $bele=$bits[0];
   if($ele eq $bele) {$skip=$.+1};
   if($skip == $.) {$writit=1};
   if($writit == 1) {print TMP "$line2\n"};
   }
  }
  close GAMB;
  print TMP "\n";
}else{
  close TMP;
  system("cat gbfile >> tmpfile");
  open(TMP,">>tmpfile");
}
}
  print TMP "\$END\n";
  close TMP;
  close TMP1;
  system("mv tmpfile FRAG$ic.inp");
#  print IN "FRAG$ic.inp\n";
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

  $ch=`awk 'NR==3,NR==3{print \$1}' ab.$ic1.$ic2.com`;
  $mult=`awk 'NR==3,NR==3{print \$2}' ab.$ic1.$ic2.com`;
  chomp $ch;
  chomp $mult;
  $mult=~ s/^\s+|\s+$//g;
  open(TMP,">tmpfile");
  if($mult == 1) {
  print TMP " \$CONTRL $prop MULT=$mult ICHARG=$ch \$END\n";
  }
  if($mult == 2) {
  print TMP " \$CONTRL $prop2 MULT=$mult ICHARG=$ch \$END\n";
  }
  close TMP;

  if ( -s "GAMextras") {
   if( $deriv == 2 && $oldderiv > 2 ) {
# this is the hessian at the start of an opt or ts
# so we only use scf for gamess
    system("grep SYSTEM GAMextras >> tmpfile");
   }else{
   system("cat GAMextras >> tmpfile");
   }
  }

  open(TMP,">>tmpfile");
  print TMP " \$DATA\n";
# get the sign of the ab nb combination
  $_=`head -1 ab.$ic1.$ic2.com`;
  s/!Isg_coeff=/\$comment!Isg_coeff=/;
  $comm=$_;
  chomp $comm;
  print TMP " $comm\n";
  print TMP "c1\n";
  my $lf=`wc -l ab.$ic1.$ic2.com`;
  chomp $lf;
  my @kf=split(/\s/,"$lf");
  $mf=$kf[0];
  $mf=$mf-3;
  system("tail -$mf ab.$ic1.$ic2.com > tfile1");
  $mf=$mf-1;
  system("head -$mf tfile1 > tfile");
  system("$EXEDIR/convertcoords4games < tfile > gamesscoords");
  open(TMP1,"<gamesscoords");
  while (<TMP1>) {
   my $line1=$_;
   chomp $line1;
   my @bits=split(/\s+/,$line1);
   $ele=$bits[0];
   print TMP "$line1\n";
if ($ans eq "N") {
  open(GAMB,"<gbfile");
  $writit=0;
  $skip=0;
  while (<GAMB>) {
   my $line2=$_;
   chomp $line2;
   if($line2 eq "") {
    $writit=0;
   }else{
   my @bits=split(/\s+/,$line2);
   $bele=$bits[0];
   if($ele eq $bele) {$skip=$.+1};
   if($skip == $.) {$writit=1};
   if($writit == 1) {print TMP "$line2\n"};
   }
  }
  close GAMB;
  print TMP "\n";
}else{
  close TMP;
  system("cat gbfile >> tmpfile");
  open(TMP,">>tmpfile");
}
} 
  print TMP "\$END\n";
  close TMP;
  close TMP1;
  system("mv tmpfile ab.$ic1.$ic2.inp");
#  print IN "ab.$ic1.$ic2.inp\n";
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
  $ch=`awk 'NR==1,NR==1{print \$1}' nb.$ic1.0.com`;
  $mult=`awk 'NR==1,NR==1{print \$2}' nb.$ic1.0.com`;
  chomp $ch;
  chomp $mult;
  $mult=~ s/^\s+|\s+$//g;
  open(TMP,">tmpfile");
  if($mult == 1) {
  print TMP " \$CONTRL $prop MULT=$mult ICHARG=$ch \$END\n";
  }
  if($mult == 2) {
  print TMP " \$CONTRL $prop2 MULT=$mult ICHARG=$ch \$END\n";
  }
  close TMP;

  if ( -s "GAMextras") {
   if( $deriv == 2 && $oldderiv > 2 ) {
# this is the hessian at the start of an opt or ts
# so we only use scf for gamess
    system("grep SYSTEM GAMextras >> tmpfile");
   }else{
   system("cat GAMextras >> tmpfile");
   }
  }

  open(TMP,">>tmpfile");
  print TMP " \$DATA\n";
  print TMP " GAMESS job\n";
  print TMP "c1\n";


  my $lf=`wc -l nb.$ic1.0.com`;
  chomp $lf;
  my @kf=split(/\s/,"$lf");
  $mf=$kf[0];
  $mf=$mf-1;
  system("tail -$mf nb.$ic1.0.com > tfile1");
  $mf=$mf-1;
  system("head -$mf tfile1 > tfile");
  system("$EXEDIR/convertcoords4games < tfile > gamesscoords");
  open(TMP1,"<gamesscoords");
  while (<TMP1>) {
   my $line1=$_;
   chomp $line1;
   my @bits=split(/\s+/,$line1);
   $ele=$bits[0];
   print TMP "$line1\n";
if ($ans eq "N") {
  open(GAMB,"<gbfile");
  $writit=0;
  $skip=0;
  while (<GAMB>) {
   my $line2=$_;
   chomp $line2;
   if($line2 eq "") {
    $writit=0;
   }else{
   my @bits=split(/\s+/,$line2);
   $bele=$bits[0];
   if($ele eq $bele) {$skip=$.+1};
   if($skip == $.) {$writit=1};
   if($writit == 1) {print TMP "$line2\n"};
   }
  }
  close GAMB;
  print TMP "\n";
}else{
  close TMP;
  system("cat gbfile >> tmpfile");
  open(TMP,">>tmpfile");
}
}
  print TMP "\$END\n";
  close TMP;
  close TMP1;
  system("echo ' \$stone' >> tmpfile");
  system("echo '  bigexp=0.0' >> tmpfile");
  system("echo '  atoms' >> tmpfile");
  system("echo ' \$end' >> tmpfile");
  system("mv tmpfile nb.$ic1.0.inp");
#  print IN "nb.$ic1.0.inp\n";
  }
close $fh;

# the nb polar files are constructed for DALTON
# as  for NWChem 
# see rundalpolar
#close IN;
}

sub SMFArunallgam {

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

$file=join("",$stem,".inp");
system("rungms $file \$PBS_NCPUS > $stem.log");
$ok=`grep 'TERMINATED NORMALLY' $stem.log | wc -l`;
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














