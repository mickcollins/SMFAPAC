#!/usr/bin/perl

$CODEDIR=`awk 'NR==1,NR==1 {print \$1}' CODENAME`;
chomp($CODEDIR);
$CODEDIR=~ s/^\s+|\s+$//g;

$EXEDIR=`awk 'NR==1,NR==1 {print \$1}' EXENAME`;
chomp($EXEDIR);
$EXEDIR=~ s/^\s+|\s+$//g;

&anal_MAC;

$junk=`rm theseforces*`;
$junk=`rm thesedipdr*`;
$junk=`rm thesecoords*`;



sub anal_MAC {

#system("clear");
#system("stty sane");

      $nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
      $deriv =`awk 'NR==2,NR==2 {print \$1}' IN_JOBTYPE`;
      $fraglvl =`awk 'NR==2,NR==2 {print \$1}' IN_FRACT`;

chomp $nchgs;
chomp $deriv;
chomp $fraglvl;

system("rm -f combined* ");

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
close $fh;

$qmprog=$store[1];
for (my $ic=0; $ic<=$n; $ic++) {
 if ("$store[$ic]" eq "Is long range dispersion accounted for?")  {
  $ie=$ic+1;
  $disp=$store[$ie];
 }
}

if ($qmprog == 1) {

# &extract_gam;

&makenbgdma_gam;

# always read polar for GAMESS
 &ReadPolar_dal;
}


if ($qmprog == 2) {

# &extract_gau;

 &makenbgdma_gau;

if ("$disp" eq "Y") {
 &ReadPolar_gau;
 } else {
   if ($nchgs eq 0) {
    &ReadPolar_gau;
   }
}
}

if ($qmprog == 3) {

# &extract_nwc;

&makenbgdma_nwc;

if ("$disp" eq "Y") {
 &ReadPolar_dal;
 } else {
   if ($nchgs eq 0) {
    &ReadPolar_dal;
   }
}
}

if ($qmprog == 4) {

# &extract_qch;

&makenbgdma_qch;

if ("$disp" eq "Y") {
 &ReadPolar_qch;
 } else {
   if ($nchgs eq 0) {
    &ReadPolar_qch;
   }
}
}

if ( ("$nchgs" eq "0") or ( $qmprog == 1 ) ) {
 system("$EXEDIR/Induct_Mac");
}

 if($disp eq "Y") {
  &extract_disp_dal;
  system("$EXEDIR/Disp");
 }



system("$EXEDIR/Derivs_smfa < FragDerivatives > rubd");
if($deriv >= 1) {
system("$EXEDIR/ABNBderivatives");
}
system("$EXEDIR/Elect_Mac > rube");
if ( ($nchgs > 1) && ( $qmprog != 1)) {
system("$EXEDIR/DoubleC_Mac > rubD");
}
system("$EXEDIR/combineallderivatives");
system("cat SMFA.out.temp >> OUT_SMFA");
system("rm SMFA.out.temp");

if($deriv == 1) {
system("echo '' >> OUT_SMFA");
system("echo ' The energy gradients are contained in the file combinedderivs' >> OUT_SMFA");
}

if($deriv == 2) {
system("$EXEDIR/DipoleDerivs < FragDipderivs > combDipDerivs");
system("rm -f IN_ISOTOPES");
system("$EXEDIR/Freq");
system("echo '' >> OUT_SMFA");
system("echo 'The normal mode frequencies are printed in file FREQUENCIES'  >> OUT_SMFA");
system("echo '' >> OUT_SMFA");
system("echo 'The normal mode eigenvectors are printed in file NORMAL_MODES' >> OUT_SMFA");
system("echo '' >> OUT_SMFA");
system("echo 'A simulated IR spectrum is printed in file SPECTRUM'  >> OUT_SMFA");
system("echo '' >> OUT_SMFA");
system("grep 'The zero-point energy' FREQUENCIES >> OUT_SMFA");
}

#&totaltimes;

system("stty sane");
#system(sleep 5);
}

sub totaltimes {
if ( -s "TIMELOADING") {
@tottime=0;
open(INT,"<TIMELOADING");
while (<INT>) {
$line=$_;
chomp($line);
@timearray=split(' ',$line);
$tottime[$timearray[0]]=$tottime[$timearray[0]] + $timearray[1];
}
close INT;
$np=$#tottime;
open(OUT,">PROCTIMES");
for ($i=1;$i <= $np;$i++) {
print OUT "$i     $tottime[$i]\n";
}
close OUT;
}
}

# superceded sub maybe not
sub make_gdma
{
        my $fh=$_[0];
        print $fh "Title \"$prefix\"\n";

#Check if we need a density keyword
        foreach $GLine (@GDMALines)
        {
           if ($GLine =~ /Density/) #Density line needs to go before File if it's there
           {
              print $fh "$GLine";
           }
        }
        print $fh "File $prefix.fchk\n\n";
        print $fh "Angstrom\n\n";
        foreach $GLine (@GDMALines)
        {
                if ($GLine =~ /Punch/) #find the checkpoint file line
                {
                        chomp $GLine;
                        print $fh " Punch $prefix.punch\n";
                }
                elsif ($GLine =~ /ADD/)#origin line if there is one
                {
                        print $fh @Names;
                }
                elsif ($GLine =~ /Density/) #Density line needs to go before File if it's there
                {
                   #Do nothing

                }
                else
                {
                        print $fh "$GLine";
                }
        }
}




sub maketotalpolar {
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
open(TMP,"<IN_JOBTYPE");
while (<TMP>) {
if ($. == 2) {$deriv="$_"};
}
chomp $deriv;

$qmprog=$store[1];

# first use the group polarizabilities

if (-s "Polarisabilities.out") {
for (my $ic=1;$ic <= 7;$ic++) {
$totpolar[$ic]=0;
}

open(TMP,"<Polarisabilities.out");
while (<TMP>) {
$line1="$_";
@line=split(/\s+/,$line1);
for (my $ic=1;$ic <= 7;$ic++) {
 $totpolar[$ic]=$totpolar[$ic]+$line[$ic];
}
}
close TMP;
open(TMP,">>OUT_SMFA");
print TMP "\n";
print TMP "The total molecular dipole polarizability is approximated\n";
print TMP "by a sum over group polarizabilities as:\n";
print TMP "\n";
print TMP "     xx         yx         yy         zx         zy         zz\n";
print TMP "\n";
printf TMP "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",$totpolar[1],$totpolar[2],$totpolar[3],$totpolar[4],$totpolar[5
],$totpolar[6];
print TMP "\n";
printf TMP "%s   %10.3f\n",'The isotropic polarizability = ',$totpolar[7];
print TMP "\n";
}

if (($qmprog == 2) || ($qmprog == 4)) {
if ($deriv == 2) {
print TMP "For GAUSSIAN and QCHem, if frequencies have been calculated,\n";
print TMP "SMFA can calculate the dipole polarizability as a sum over the fragment polarizabilities.\n";
print TMP "This is a more accurate approximation.\n";
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
$totpolar[$ic]=0;
}
for (my $ic=1;$ic <= $nfrgs;$ic++) {
if($qmprog == 2) {&getpolar_gau("FRAG$ic.log")};
if($qmprog == 4) {&getpolar_qch("FRAG$ic.log")};
if($qmprog == 2) {
for (my $ie=1;$ie <= 6;$ie++) {
 $totpolar[$ie]=$totpolar[$ie]+$sign[$ic]*$Polar[$ie];
}
$totpolar[7]=$totpolar[7]+$sign[$ic]*($Polar[1]+$Polar[3]+$Polar[6]);
}
if($qmprog == 4) {
for (my $ie=1;$ie <= 9;$ie++) {
 $totpolar[$ie]=$totpolar[$ie]+$sign[$ic]*$Polar[$ie];
}
$totpolar[10]=$totpolar[10]+$sign[$ic]*($Polar[1]+$Polar[5]+$Polar[9]);
}
}
if($qmprog == 2) {
$totpolar[7]=$totpolar[7]/3;
}
if($qmprog == 4) {
$totpolar[10]=$totpolar[10]/3;
}
print TMP "\n";
print TMP "     xx         yx         yy         zx         zy         zz\n";
if($qmprog == 2) {
printf TMP "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",$totpolar[1],$totpolar[2],$totpolar[3],$totpolar[4],$totpolar[5],$totpolar[6];
print TMP "\n";
printf TMP "%s   %10.3f\n",'The isotropic polarizability = ',$totpolar[7];
}
if($qmprog == 4) {
# reverse an apparent sign error in QChem
for (my $ic=1;$ic <= 10;$ic++) {
$totpolar[$ic]=-$totpolar[$ic];
}

printf TMP "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",$totpolar[1],$totpolar[4],$totpolar[5],$totpolar[7],$totpolar[8],$totpolar[9];
print TMP "\n";
printf TMP "%s   %10.3f\n",'The isotropic polarizability = ',$totpolar[10];
print TMP "\n";
}
}else{
print TMP "A frequency calculation was not the last procedure, aborting\n";
}
}else{
print TMP "only GAUSSIAN or QChem are allowed here, aborting\n";
}

close TMP;
print "The polarizability has been evaluated and the result has been appended to OUT_SMFA\n";
system("sleep 4");

}


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
$deriv=<STDIN>;
chomp $deriv;

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

print "\n";
system("rm GAMbasis");
print "GAMESS allows a different basis set for each chemical element\n";
print "Do you want all elements to have the same basis (Y or N) ? ";
my $ans=<STDIN>;
chomp $ans;
$ans=uc($ans);
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

# get CONTRL group commands
open(TMP,"<GAMcommands");
$ic=0;
while(<TMP>) {
my $line="$_";
chomp $line;
$ic=$ic+1;
$prop=join(" ",$prop,$line);
}

  $ans=`awk 'NR==1,NR==1 {print \$1}' GAMbasis`;
  chomp $ans;
  my $lf=`wc -l GAMbasis`;
  chomp $lf;
  my @kf=split(/\s/,"$lf");
  $mf=$kf[0];
  $mf=$mf-1;
  system("tail -$mf GAMbasis > gbfile");

open(IN,">INLIST");

# start with the main fragments
  my $nfrgs = `ls COORD* | wc -l`;
  chomp $nfrgs;
 for (my $ic=1; $ic<=$nfrgs; $ic++)
  {
  $ch=`awk 'NR==1,NR==1{print \$1}' COORD$ic`;
  $mult=`awk 'NR==1,NR==1{print \$3}' COORD$ic`;
  chomp $ch;
  chomp $mult;
  open(TMP,">tmpfile");
  print TMP " \$CONTRL $prop MULT=$mult ICHARG=$ch \$END\n";
  close TMP;
  if ( -s "GAMextras") {system("cat GAMextras >> tmpfile")};
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
  system("$EXEDIR/convertcoords4games < tfile");
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
  print IN "FRAG$ic.inp\n";
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
  open(TMP,">tmpfile");
  print TMP " \$CONTRL $prop MULT=$mult ICHARG=$ch \$END\n";
  close TMP;
  if ( -s "GAMextras") {system("cat GAMextras >> tmpfile")};
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
  system("$EXEDIR/convertcoords4games < tfile");
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
  print IN "ab.$ic1.$ic2.inp\n";
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
  open(TMP,">tmpfile");
  print TMP " \$CONTRL $prop MULT=$mult ICHARG=$ch \$END\n";
  close TMP;
  if ( -s "GAMextras") {system("cat GAMextras >> tmpfile")};
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
  system("$EXEDIR/convertcoords4games < tfile");
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
  print IN "nb.$ic1.0.inp\n";
  }
close $fh;

# the nb polar files are constructed for DALTON
# as  for NWChem 
# see rundalpolar
close IN;
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
($stem)=split(".inp",$file);

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
}














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
$deriv=<STDIN>;
chomp $deriv;
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
print "Does this method account for electronic correlation leading to dispersion?\n";
print "(See the user's manual). Do you want to account for dispersion at long range?\n";
print "Answer Y or N\n";
$disp=<STDIN>;
chomp $disp;
$disp=uc($disp);

print "\n";
print "GAUSSIAN allows a different basis set for each chemical element.\n";
print "Do you want all elements to have the same basis (Y or N) ? \n";
$ans=<STDIN>;
chomp $ans;
$ans=uc($ans);
if ($ans eq "Y") {
print "Enter the basis for all atoms\n";
$basis=<STDIN>;
chomp $basis;
}else{
print "Enter the basis (check availability in the GAUSSIAN manual)\n";
print "for the following elements\n";
open (TMP,"<labels0");
system("rm GAUbasis");
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

# Gaussian job header
       open G09IN, ">IN_G09";
#       # Gaussian job Model - method/basis
#             print G09IN "\%nproc=$nproc0\n";
if ($ans eq "Y") {
             print G09IN "\#p $method/$basis nosymm density=current $prop $conv $opts\n\n";
}else{
             print G09IN "\#p $method/gen nosymm density=current $prop  $conv $opts\n\n";
}
             close G09IN;
       open G09IN, ">IN_NPA";
$popnpa="pop=NPA";
# header for NPA charges
#             print G09IN "\%nproc=$nproc0\n";
if ($ans eq "Y") {
print G09IN "\#p $method/$basis charge nosymm $popnpa\n\n";
}else{
             print G09IN "\#p $method/gen charge nosymm $popnpa\n\n";
}
             close G09IN;

       open G09IN, ">IN_POLAR";
#             print G09IN "\%nproc=$nproc0\n";
if ($ans eq "Y") {
print G09IN "\#p $method/$basis polar nosymm\n\n";
}else{
print G09IN "\#p $method/gen polar nosymm\n\n";
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
   $m1 = ltrim($n1);
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

open(IN,">INLIST");

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
print IN "FRAG$ic.com\n";
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
 print IN "ab.$ic1.$ic2.com\n";
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
 print IN "nb.$ic1.0.com\n";
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
 print IN "nb.$ic1.0-polar.com\n";
}
close $fh;
close IN;
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
$deriv=<STDIN>;
chomp $deriv;
#if($deriv==0){$prop=" "};
#if($deriv==1){$prop="FORCE"};
#if($deriv==2){$prop="FREQ"};
#if($deriv==3){$prop="OPT"};
#if($deriv==4){$prop="OPT=TS"};

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
#print "Enter the number of processors for each calculation:\n ";
#$nproc0=<STDIN>;
#chomp $nproc0;

print "\n";

print "NWChem allows a different basis set for each element\n";
print "Do you want all elements to have the same basis (Y or N) ? ";
$ans=<STDIN>;
chomp $ans;
$ans=uc($ans);
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
system("rm NWCbasis");
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

open(IN,">INLIST");

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
  if( $mult == 2 ) {
   s/SINGLET/$mform/;
  }  
  print TMP "$_";
 }
  close TMP;
  system("mv tmpfile FRAG$ic.nw");
  print IN "FRAG$ic.nw\n";
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
  if( $mult == 2 ) {
   s/SINGLET/$mform/;
  }
  print TMP "$_";
 }
  close TMP;
  system("mv tmpfile ab.$ic1.$ic2.nw");
  print IN "ab.$ic1.$ic2.nw\n";
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
  if( $mult == 2 ) {
   s/SINGLET/$mform/;
  }
  print TMP "$_";
 }
  print TMP "ESP\n";
  print TMP "recalculate\n";
  print TMP "END\n";
  print TMP "TASK ESP";
  close TMP;
  system("mv tmpfile nb.$ic1.0.nw");
  print IN "nb.$ic1.0.nw\n";
}
close $fh;

# the nb polar files are constructed for DALTON
# as NWChem does not have polarizabilities
# see rundalpolar
close IN;
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
  print TMP " direct";
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
  print TMP "TASK ENERGY SCF\n";
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
($stem)=split(".nw",$file);

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
system ("ln -s IN_DISP $stem.dal");
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
}















sub readqchparams {

if ( -s "xyzFILENAME" ) {
$file=`awk 'NR==1,NR==1 {print \$1}' xyzFILENAME`;
}else{
print "You must first enter the filename for the Cartesian coordinates\n";
system ("sleep 5");
return 1;
}
system("$EXEDIR/distinctlabels < $file ");
system("mv labels labels0");

print "                       Variables needed by SMFA and QChem\n";
print "\n";


print "Enter values for the \$rem QChem input\n";
print "The JOBTYPE can be SP, FORCE, FREQ, OPT, TS or SCAN\n";
print "JOBTYPE ";
$jobtype=<STDIN>;
chomp $jobtype;
$jobtype=uc($jobtype);
print "METHOD ";
$method=<STDIN>;
chomp $method;
$method=uc($method);
if($method eq "SCF") {$method="HF"};
print "Enter any additional values for the \$rem QChem input\n";
print "Do not enter a BASIS value\n";
print "Begin now, and end with a RETURN\n";

my $ic=0;
my $again="Y";
while ($again ne "N") {
my $inp=<STDIN>;
if ($inp eq "\n") {
 $again="N";
}else{
chomp $inp;
$ic=$ic+1;
$qchcommands[$ic]=uc($inp);
}
}

print "Does this METHOD account for electronic correlation leading to dispersion?\n";
print "(See the user's manual). If YES, then do you want to account for dispersion at long range?\n";
print "Answer Y or N\n";
$disp=<STDIN>;
chomp $disp;
$disp=uc($disp);
if ("$jobtype" eq "SP") {$deriv=0}
if ("$jobtype" eq "FORCE") {$deriv=1}
if ("$jobtype" eq "FREQ") {$deriv=2}
if ("$jobtype" eq "OPT") {$deriv=3}
if ("$jobtype" eq "TS") {$deriv=4}
if ("$jobtype" eq "SCAN") {$deriv=5}
print "\n";

print "QChem allows a different basis set for each element\n";
print "Do you want all elements to have the same basis (Y or N) ? ";
my $ans=<STDIN>;
chomp $ans;
$ans=uc($ans);
if ($ans eq "Y") {
print "Enter the basis for all atoms\n";
$basis=<STDIN>;
chomp $basis;
}else{
print "Enter the basis (check availability in the QChem library)\n";
print "for the following elements\n";
open (TMP,"<labels0");
system("rm QCHbasis");
open (TMP1,">QCHbasis");
print TMP1 "\$basis\n";
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
print TMP1 "\$end\n";
close TMP;
close TMP1;
}

&solventcharges;


if($deriv==3) {&optinput};
if($deriv==4) {&optinput};
if($deriv==4) {&TSinput};
if($deriv==5) {&scaninput};

#print "\n";
#print "Enter the number of processors for each calculation:\n ";
#$nproc0=<STDIN>;
#chomp $nproc0;

open TMP,">IN_QCH";
print TMP "\$rem\n";
if($deriv==5) {
$jobtype="OPT"
} 
if($deriv==4) {
$jobtype="OPT"
} 
print TMP "JOBTYPE $jobtype\n";
if($deriv >= 2) {print TMP "IPRINT 10000000\n"};
print TMP "METHOD $method\n";
if ($ic > 0) {
for (my $ie=1;$ie <= $ic;$ie++) {
 print TMP "$qchcommands[$ie]\n";
}
}
if ($ans eq "Y") {
print TMP "BASIS $basis\n";
}else{
print TMP "BASIS General\n";
print TMP "PURECART 2\n";
}
print TMP "GUI 2\n";
print TMP "SYM_IGNORE TRUE\n";
print TMP "\$end\n";
close TMP;
if ($ans eq "N") {
system("cat QCHbasis >> IN_QCH");
}

open TMP,">IN_QCHNPA";
print TMP "\$rem\n";
print TMP "JOBTYPE SP\n";
print TMP "METHOD $method\n";
if ($ic > 0) {
for (my $ie=1;$ie <= $ic;$ie++) {
 print TMP "$qchcommands[$ie]\n";
}
}
if ($ans eq "Y") {
print TMP "BASIS $basis\n";
}else{
print TMP "BASIS General\n";
print TMP "PURECART 2\n";
}
print TMP "GUI 2\n";
print TMP "SYM_IGNORE TRUE\n";
print TMP "NBO ON\n";
print TMP "\$end\n";
close TMP;
if ($ans eq "N") {
system("cat QCHbasis >> IN_QCHNPA");
}


open TMP,">IN_QCHPOLAR";
print TMP "\$rem\n";
print TMP "JOBTYPE SP\n";
print TMP "METHOD $method\n";
if ($ic > 0) {
for (my $ie=1;$ie <= $ic;$ie++) {
 print TMP "$qchcommands[$ie]\n";
}
}
if ($ans eq "Y") {
print TMP "BASIS $basis\n";
}else{
print TMP "BASIS General\n";
print TMP "PURECART 2\n";
}
print TMP "GUI 2\n";
print TMP "SYM_IGNORE TRUE\n";
print TMP "MOPROP 2\n";
print TMP "\$end\n";
close TMP;
if ($ans eq "N") {
system("cat QCHbasis >> IN_QCHPOLAR");
}


      open TMP, ">ABSETTINGS";
       print TMP "The quantum chemistry program package is\n";
       print TMP "$package\n";
       print TMP "The QChem input variables are\n";
       close TMP;
       system("cat IN_QCH >> ABSETTINGS");
      open TMP, ">>ABSETTINGS";
       print "The basis set\(s\)\n";
       close TMP;
       system("cat QCHbasis >> ABSETTINGS");
      open TMP, ">>ABSETTINGS";
       print TMP "The job type is energy (0), force (1), hessian (2), Opt (3), TS (4) or scan (5)\n";
       print TMP "$deriv\n";
       print TMP "Is long range dispersion accounted for?\n";
       print TMP "$disp\n";
#       print TMP "The number of processors for individual calculations is\n";
#       print TMP "$nproc0\n";
      close TMP;

# make the IN_GDMA file
open(TMP,">IN_GDMA");
print TMP "Density SCF\n";
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

# make the IN_JOBTYPE file
open(TMP,">IN_JOBTYPE");
print TMP "Enter 0, 1, 2, 3, 4, 5 for energy, gradient, hessian, opt, TS or scan\n";
print TMP "$deriv\n";
close TMP;

}

sub extract_qch {

# get the energy, force and hessian from all re;evant ab initio jobs
# run by QChem

system("clear");
system("stty echo");

# some data need in getcoords_qch
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

&getcoords_qch("FRAG$ic.log");
system("cat thesecoords >> FragDerivatives");

&getenergy_qch("FRAG$ic.log");
open(FR,">>FragDerivatives");
print FR "$TotEn\n";
close FR;
  if ($deriv >= 1) {
   &getforces_qch("FRAG$ic.log");
   system("cat theseforces >> FragDerivatives");
  }
  if ($deriv >= 2) {
   &gethessian_qch("FRAG$ic.com.fchk");
   system("cat hessian >> FragDerivatives");
   &getdipdr_qch("FRAG$ic.log");
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
 &getenergy_qch("ab.$ic1.$ic2.log");
 my $TotE=$TotEn;
 &getenergy_qch("nb.$ic1.0.log");
 $TotE=$TotE - $TotEn;
 &getenergy_qch("nb.$ic2.0.log");
 $TotE=$TotE - $TotEn;
# $absign=`awk '/Isg_coeff/ {print \$2}' ab.$ic1.$ic2.com`;
 system("awk /coeff/ ab.$ic1.$ic2.com > rub");
 $absign=`awk 'NR==1,NR==1 {print \$2}' rub`;
 chomp $absign;
 $AllTot = $AllTot + $absign * $TotE;

 if ($deriv >= 1) {
  print TMP1 "$ic1","  ","$ic2\n";
  print TMP1 "$absign\n";
  &getforces_qch("ab.$ic1.$ic2.log");
  print TMP1 "$numberofatoms\n";
  my $nat1=$numberofatoms;
  close TMP1;
  system("cat theseforces >> abforces");
  open(TMP1,">>abforces");
  print TMP1 "$ic1\n";
  &getforces_qch("nb.$ic1.0.log");
  print TMP1 "$numberofatoms\n";
  my $nat2=$numberofatoms;
  close TMP1;
  system("cat theseforces >> abforces");
  open(TMP1,">>abforces");
  print TMP1 "$ic2\n";
  &getforces_qch("nb.$ic2.0.log");
    print TMP1 "$numberofatoms\n";
    my $nat3=$numberofatoms;
  close TMP1;
  system("cat theseforces >> abforces");
  open(TMP1,">>abforces");


 if ($deriv >= 2) {
  print TMP2 "$ic1","  ","$ic2\n";
  print TMP2 "$absign\n";
  &gethessian_qch("ab.$ic1.$ic2.com.fchk");
  print TMP2 "$nat1\n";
  close TMP2;
  system("cat hessian >> abhessians");
  open(TMP2,">>abhessians");
  print TMP2 "$ic1\n";
  &gethessian_qch("nb.$ic1.0.com.fchk");
  print TMP2 "$nat2\n";
  close TMP2;
  system("cat hessian >> abhessians");
  open(TMP2,">>abhessians");
  print TMP2 "$ic2\n";
  &gethessian_qch("nb.$ic2.0.com.fchk");
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

sub getcoords_qch {
 open(FRAGJOB,"<$_[0]");
 open(TMP,">thesecoords");
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
open(TMP,">thesedipdr");
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
      $entc = "$totline[5]";
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
open (TMP,">tempfile");
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
open (TMP,">theseforces");
my $line1=`awk 'NR==1,NR==1 {print \$0}' tempfile`; 
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
 my $line1=`awk 'NR=='$newline',NR=='$newline' {print \$0}' tempfile`;
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
 my $line1=`awk 'NR=='$newline',NR=='$newline' {print \$0}' tempfile`;
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
system("rm tempfile");

}

sub gethessian_qch {

my $filename="$_[0]";
# filename is an fchk file
my $numberofatoms=`awk '/Number of atoms/ {print \$5}' $filename`;
open (FH,"<$filename");
open (TMP,">hessian");
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


sub makenbgdma_qch {
open(GDMA,"IN_GDMA") or die "Unable to open IN_GDMA";
@GDMALines = <GDMA>;

$pwd=`pwd`;
chomp $pwd;

@All_jobs= `ls -lah nb*.0.com.fchk --format=single-column`;
         $Number_of_jobs=@All_jobs;
         foreach $job (@All_jobs){
            chomp $job;
            $prefix=$job;
            $prefix =~ s/\.[^.]*$//;

            @bits=split(/\.com/,"$prefix");
            $prefix=$bits[0];
system("mv $prefix.com.fchk $prefix.fchk");

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
          }
}


sub ReadPolar_qch  {
system("rm Polarisabilities.out");
open(my $OUT, ">Polarisabilities.out");

my $nfrgs = `ls nb.*.0-polar.log | wc -l`;
for (my $ic=1; $ic<=$nfrgs; $ic++) {
$job="nb.$ic.0-polar.log";
&getpolar_qch("$job");
$iso=($Polar[1] + $Polar[5] + $Polar[9])/3;
print $OUT "$ic\t$Polar[1]\t$Polar[4]\t$Polar[5]\t$Polar[7]\t$Polar[8]\t$Polar[9]\t$iso\n";

}
}

sub getpolar_qch {
 my $filename="$_[0]";
 open(FRAGJOB,"<$_[0]");
 $skip=0;
 $skip2=0;
 $writit=0;
 $nvar=1;
 while (<FRAGJOB>) {
 if (/Polarizability \(a.u.\)/) {
   $skip=$. + 15;
   $skip2=$. + 18;}
 if (/Polarizability Matrix \(a.u.\)/) {
   $skip=$. + 2;
   $skip2=$. + 5;}
if ($. == $skip) {$writit=1};
if ($. == $skip2) {$writit=0};
      if( $writit == 1) {
       @line=split(/\s+/,"$_");
       $l3=$#line;
       $l2=$l3-1;
       $l1=$l2-1;
       $Polar[$nvar]=$line[$l1];
       $nvar=$nvar+1;
       $Polar[$nvar]=$line[$l2];
       $nvar=$nvar+1;
       $Polar[$nvar]=$line[$l3];
       $nvar=$nvar+1;

      }
 }
 close FRAGJOB;

}


sub SMFArunqch {

$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
$n++;
}
close $fh;
$n=$n-1;
#$nproc0=$store[$n];

# temporary 250118
$nproc0=1;

for (my $ic=0; $ic<=$n; $ic++) {
 if ("$store[$ic]" eq "Is long range dispersion accounted for?")  {
  $ie=$ic+1;
  $disp=$store[$ie];
 }
}

  my $nfrgs = `ls COORD* | wc -l`;
chomp $nfrgs;

 for (my $ic=1; $ic<=$nfrgs; $ic++)
  {
  system("qchem -nt $nproc0 FRAG$ic.com FRAG$ic.log");
 $ok=`grep 'Have a nice day' FRAG$ic.log | wc -l`;
 chomp($ok);
 $ok=~ s/^\s+|\s+$//g;
 if ( $ok == 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation failed for FRAG$ic\n";
  close TMP;
  exit(0);
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
 system("qchem -nt $nproc0 ab.$ic1.$ic2.com ab.$ic1.$ic2.log");
 $ok=`grep 'Have a nice day' ab.$ic1.$ic2.log | wc -l`;
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
 system("qchem -nt $nproc0 nb.$ic1.0.com nb.$ic1.0.log");
 $ok=`grep 'Have a nice day' nb.$ic1.0.log | wc -l`;
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
 system("qchem -nt $nproc0 nb.$ic1.0-polar.com nb.$ic1.0-polar.log");
 $ok=`grep 'Have a nice day' nb.$ic1.0-polar.log | wc -l`;
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
 system("qchem -nt $nproc0 nb.$ic1.0-polar.com nb.$ic1.0-polar.log");
 $ok=`grep 'Have a nice day' nb.$ic1.0-polar.log | wc -l`;
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



sub SMFAqchinputs_MAC {

    $nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;

open(IN,">INLIST");

# start with the main fragments
  my $nfrgs = `ls COORD* | wc -l`;
 chomp $nfrgs;

 for (my $ic=1; $ic<=$nfrgs; $ic++)
  {
  my $my_coord = "COORD${ic}";
  system("cp IN_QCH FRAG$ic.com");
  system("echo '\$molecule' >> FRAG$ic.com");
  system("cat $my_coord >> FRAG$ic.com");
  system("echo '\$end' >> FRAG$ic.com");
  if ( -s "chFRAG.$ic")  {
   system("echo '\$external_charges' >> FRAG$ic.com");
   open (CHID,"chFRAG.$ic");
   foreach my $n1 (<CHID>) {
   chomp $n1;
   $m1 = &ltrim($n1);
   system("cat charge.$m1.grp >> FRAG$ic.com");
   }
   close CHID;
   system("echo '\$end' >> FRAG$ic.com");
   }
# if deriv=2 (hessian) then we need to add a first job
# for the gradient
if($deriv == 2) {
 system("cp FRAG$ic.com tmpfile1");
 system("sed 's/JOBTYPE FREQ/JOBTYPE FORCE/' tmpfile1 > tmpfile");
 system("echo '\@\@\@' >> tmpfile");
 system("echo '' >> tmpfile");
 system("cat FRAG$ic.com >> tmpfile");
 system("mv tmpfile FRAG$ic.com");
 print IN "FRAG$ic.com\n";
}

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
 system("cp IN_QCH tmpfile");
 system("echo '\$molecule' >> tmpfile");
# get the sign of the ab nb combination
 $_=`head -1 ab.$ic1.$ic2.com`;
 s/!Isg_coeff=/\$comment!Isg_coeff=/;
 my $comm=$_;
 my $lf=`wc -l ab.$ic1.$ic2.com`;
 chomp $lf;
 my @kf=split(/\s/,"$lf");
 $mf=$kf[0];
 $mf=$mf-1;
 system("head -$mf ab.$ic1.$ic2.com > tfile");
 $mf=$mf-2;
 system("tail -$mf tfile  >> tmpfile");
 system("echo '\$end' >> tmpfile");
 system("echo '$comm' >> tmpfile");
 if($nchgs >= 1){
  system("echo '\$external_charges' >> tmpfile");
  open (CHID,"chab.$ic1.$ic2");
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $m1 = ltrim($n1);
   system("cat charge.$m1.grp >> tmpfile");
   }
  close CHID;
  system("echo '\$end' >> tmpfile");
  }
 system("mv tmpfile ab.$ic1.$ic2.com");
if($deriv == 2) {
 system("cp ab.$ic1.$ic2.com tmpfile1");
 system("sed 's/JOBTYPE FREQ/JOBTYPE FORCE/' tmpfile1 > tmpfile");
 system("echo '\@\@\@' >> tmpfile");
 system("echo '' >> tmpfile");
 system("cat ab.$ic1.$ic2.com >> tmpfile");
 system("mv tmpfile ab.$ic1.$ic2.com");
}
 print IN "ab.$ic1.$ic2.com\n";
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
 system("cp IN_QCH tmpfile");
 system("echo '\$molecule' >> tmpfile");
 my $lf=`wc -l nb.$ic1.0.com`;
 chomp $lf;
 my @kf=split(/\s/,"$lf");
 $mf=$kf[0];
 $mf=$mf-1;
 system("head -$mf nb.$ic1.0.com >> tmpfile");
 system("echo '\$end' >> tmpfile");
 if($nchgs >= 1){
  system("echo '\$external_charges' >> tmpfile");
  open (CHID,"chnb.$ic1");
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $m1 = ltrim($n1);
   system("cat charge.$m1.grp >> tmpfile");
   }
  close CHID;
  system("echo '\$end' >> tmpfile");
  }
 system("mv tmpfile nb.$ic1.0.com");
if($deriv == 2) {
 system("cp nb.$ic1.0.com tmpfile1");
 system("sed 's/JOBTYPE FREQ/JOBTYPE FORCE/' tmpfile1 > tmpfile");
 system("echo '\@\@\@' >> tmpfile");
 system("echo '' >> tmpfile");
 system("cat nb.$ic1.0.com >> tmpfile");
 system("mv tmpfile nb.$ic1.0.com");
}
 print IN "nb.$ic1.0.com\n";
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
 system("cp IN_QCHPOLAR tmpfile");
 system("echo '\$molecule' >> tmpfile");
 my $lf=`wc -l nb.$ic1.0-polar.com`;
 chomp $lf;
 my @kf=split(/\s/,"$lf");
 $mf=$kf[0];
 $mf=$mf-1;
 system("head -$mf nb.$ic1.0-polar.com >> tmpfile");
 system("echo '\$end' >> tmpfile");
 system("mv tmpfile nb.$ic1.0-polar.com");
 print IN "nb.$ic1.0-polar.com\n";
}
close $fh;
close IN;
# end of subroutine
}


sub Lev0_chg_MAC_iter_qch {

system("clear");

$nproc0=1;

$iter=3;
$count=1;
while ($count <= $iter) {

&SMFAqchnpa_MAC;
&extractch_qch;

      $count++;
}

      return 1;
}

sub extractch_qch {
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

sub SMFAqchnpa_MAC {

&SMFAmkchinputs_qch;
&runchseq_qch;

}

sub SMFAmkchinputs_qch {
 my @jobs = `ls charge.*.coord`;
 foreach my $chg_coord (@jobs) {
  @fields=split(/\./,$chg_coord);
  $ic=$fields[1];
  system("cp IN_QCHNPA charge.$ic.com");
  system("echo '\$molecule' >> charge.$ic.com");
  system("cat charge.$ic.coord >> charge.$ic.com");
  system("echo '\$end' >> charge.$ic.com");
  if ( -s "chL0.$ic")  {
   system("echo '\$external_charges' >> charge.$ic.com");
   open (CHID,"chL0.$ic");
   foreach my $n1 (<CHID>) {
   chomp $n1;
   $m1 = &ltrim($n1);
   system("cat charge.$m1.grp >> charge.$ic.com");
   }
   close CHID;
   system("echo '\$end' >> charge.$ic.com");
  }
}
# end sub
}

sub  runchseq_qch {

 my @jobs = `ls charge.*.coord`;
 foreach my $chg_coord (@jobs) {
  @fields=split(/\./,$chg_coord);
  $ic=$fields[1];
  system("qchem -nt $nproc0 charge.$ic.com charge.$ic.log");
 $ok=`grep 'Have a nice day' charge.$ic.log | wc -l`;
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

sub SMFArunallqch {

$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
$n++;
}
close $fh;
$n=$n-1;
$nproc0=$store[$n];

for (my $ic=0; $ic<=$n; $ic++) {
 if ("$store[$ic]" eq "Is long range dispersion accounted for?")  {
  $ie=$ic+1;
  $disp=$store[$ie];
 }
}

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
($stem)=split(".com",$file);

system("qchem -nt $nproc0 $file $stem.log");
$ok=`grep 'Have a nice day' $stem.log | wc -l`;
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
system ("ln -s IN_DISP $stem.dal");
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
}





sub extract_disp_dal  {
 system("rm PAA.out");
 open($PAA_output,">PAA.out");

 my $ngroups = `ls grp.*.mol | wc -l`;
  for(my $group=1;$group<=$ngroups;$group++)
  {
   my $filebase="grp.$group";
   $disp_out="$filebase.disp";
   $dalton_input="$filebase.mol"; 

      open($PAA_In, ">$disp_out");

      open($DAL_IN,"$dalton_input") or die "Unable to open $DAL_IN";
      @dalton_in = <$DAL_IN>;
      while (my $DLine = shift (@dalton_in))
      {
         if ($DLine =~ /^(\w+\s+)((-)*\d+\.\d+\s+)+/) #match  line that looks like coordinates
         {
            $num_atoms+=1;
            $mol_coords[$num_atoms]=$DLine;
         }
      }
      close $DAL_IN ;
      print $PAA_In "$num_atoms\n";
     for (my $ii=1; $ii<=$num_atoms; $ii++)
     {
      print $PAA_In $mol_coords[$ii];
     }
      undef(@mol_coords);

      @dalton_output=`grep -A 17 -B 1 "GRIDSQ" $filebase.out`;

      foreach $line (@dalton_output)
      {
         if ($line !~ /(GRIDSQ|--)/) # the "--" is output by grep wen there's a blank line
              {
                      print $PAA_In "$line";
              }
      }
      #Now get the PAA value and print it to a file
      @PAA_guff=`$EXEDIR/PAA <$filebase.disp`;
      my $last=$num_atoms+3;
      print $PAA_output "$group", "   ",  "$PAA_guff[$last]";
      $num_atoms=0;
  }
}

sub ReadPolar_dal {
system("rm Polarisabilities.out");
open(my $OUT, ">Polarisabilities.out");

my $nfrgs = `ls nb.*.0-polar.out | wc -l`;
chomp $nfrgs;
for (my $ic=1; $ic<=$nfrgs; $ic++) {
$job="nb.$ic.0-polar.out";
&getpolar_dal("$job");
$iso=($Polar[1] + $Polar[5] + $Polar[9])/3;
print $OUT "$ic\t$Polar[1]\t$Polar[4]\t$Polar[5]\t$Polar[7]\t$Polar[8]\t$Polar[9]\t$iso\n";

}
close $OUT;
}

sub getpolar_dal  {
 $file="<$_[0]";
 $ans=`grep 'FREQUENCY INDEPENDENT SECOND ORDER PROPERTIES' $file | wc -l`;
 chomp($ans);
 if ( $ans == 1) {
  $Polar[1]=`awk '/XDIPLEN  ; XDIPLEN/ {print \$8}' $file`;
  $Polar[2]=`awk '/XDIPLEN  ; YDIPLEN/ {print \$8}' $file`;
  $Polar[3]=`awk '/XDIPLEN  ; ZDIPLEN/ {print \$8}' $file`;
  $Polar[4]=$Polar[2];
  $Polar[5]=`awk '/YDIPLEN  ; YDIPLEN/ {print \$8}' $file`;
  $Polar[6]=`awk '/YDIPLEN  ; ZDIPLEN/ {print \$8}' $file`;
  $Polar[7]=$Polar[3];
  $Polar[8]=$Polar[6];
  $Polar[9]=`awk '/ZDIPLEN  ; ZDIPLEN/ {print \$8}' $file`;
  for ($i=1;$i <= 9;$i++) {
   chomp($Polar[$i]);
  }
  return;
 }

 open(FRAGJOB,"<$_[0]");
 $skip=0;
 $skip2=0;
 $writit=0;
 $nvar=1;
while (<FRAGJOB>)  {
 if(/Static polarizabilities \(au\)/) {
   $skip=$. + 5;
   $skip2=$. + 8;
}
if ($. == $skip) {$writit=1};
if ($. == $skip2) {$writit=0};
      if( $writit == 1) {
       $line1="$_";
       chomp $line1;
       @line=split(/\s+/,"$line1");
       $Polar[$nvar]=$line[2];
       $nvar=$nvar+1;
       $Polar[$nvar]=$line[3];
       $nvar=$nvar+1;
       $Polar[$nvar]=$line[4];
       $nvar=$nvar+1;
      }
 }
 close FRAGJOB;

}


sub dodaltonpolar {
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
close $fh;

# find which package
# and what electronic structure method
if($store[1] == 3) {
#its nwc
$here=`grep "TASK SCF" "ABSETTINGS" | wc -l`;
chomp $here;
if($here > 0) {
 $method="HF";
}else{
 $method="MP2";
}
}

if($store[1] == 1) {
#its gam
$here=`grep "HF" "ABSETTINGS" | wc -l`;
chomp $here;
if($here > 0) {
 $method="HF";
} 
$here=`grep "MPLEVL" "ABSETTINGS" | wc -l`;
chomp $here;
if($here > 0) {
 $method="MP2";
}
}


for (my $ic=0; $ic<=$n; $ic++) {
 if ("$store[$ic]" eq "Is long range dispersion accounted for?")  {
  $ie=$ic+1;
  $disp=$store[$ie];
 }
}

# make the input files that DALTON needs
&daltonpolar;

$nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;
if ($store[1] == 1) {$nchgs=0};

if ("$disp" eq "Y") {
&SMFAmkdalpolar;
#&SMFArundalpolar;
 } else {
   if ($nchgs eq 0) {
&SMFAmkdalpolar;
#&SMFArundalpolar;
   }
}

}


sub SMFAmkdalpolar {
system("ls -1 nb.*.0-polar.com > polarfiles");
$filename='polarfiles';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";

open(IN,">>INLISTDAL");

while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 system("$EXEDIR/convert2dalton < nb.$ic1.0-polar.com > nb.$ic1.0-polar.mol");
 print IN "nb.$ic1.0-polar.mol\n";
 $mult=`awk 'NR==1,NR==1 {print \$2}' nb.$ic1.0-polar.com`;
 chomp $mult;
 if($method eq "HF") {
  if($mult == 1) {
   system("cp IN_DAL_POLAR_HF_S nb.$ic1.0-polar.dal");
   }else{
   system("cp IN_DAL_POLAR_HF_D nb.$ic1.0-polar.dal");
   }
  }else{
    if($mult == 1) {
   system("cp IN_DAL_POLAR_MP2_S nb.$ic1.0-polar.dal");
   }else{
#   system("cp IN_DAL_POLAR_MP2_D nb.$ic1.0-polar.dal");
#   # approximation to stop dalton crashing with MP2 doublets
   system("cp IN_DAL_POLAR_HF_D nb.$ic1.0-polar.dal");
   }
  }
}
close $fh;
close IN;
}


sub SMFArundalpolar {

 system("ls -1 nb*.0-polar.mol > nbfilesp");
$filename='nbfilesp';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 system("dalton -N 1 -o nb.$ic1.0-polar.out nb.$ic1.0-polar");
 }
close $fh;
}

sub daltonpolar {

open(TMP,">IN_DAL_POLAR_HF_D");
print TMP "**DALTON INPUT\n";
print TMP ".RUN PROPERTIES\n";
print TMP "**WAVE FUNCTIONS\n";
print TMP ".HF\n";
print TMP "*CONFIGURATION INPUT\n";
print TMP ".SPIN MULTIPLICITY\n";
print TMP "  2\n";
print TMP "**PROPERTIES\n";
print TMP ".POLARI\n";
print TMP "**END OF DALTON INPUT\n";
close TMP;

open(TMP,">IN_DAL_POLAR_HF_S");
print TMP "**DALTON INPUT\n";
print TMP ".RUN PROPERTIES\n";
print TMP "**WAVE FUNCTIONS\n";
print TMP ".HF\n";
print TMP "**PROPERTIES\n";
print TMP ".POLARI\n";
print TMP "**END OF DALTON INPUT\n";
close TMP;

open(TMP,">IN_DAL_POLAR_MP2_D");
print TMP "**DALTON INPUT\n";
print TMP ".RUN RESPONSE\n";
print TMP "**WAVE FUNCTIONS\n";
print TMP ".HF\n";
print TMP ".MP2\n";
print TMP "*CONFIGURATION INPUT\n";
print TMP ".SPIN MULTIPLICITY\n";
print TMP "  2\n";
print TMP "**RESPONSE\n";
print TMP ".SOPPA\n";
print TMP "*LINEAR\n";
print TMP ".DIPLEN\n";
print TMP ".FREQUENCIES\n";
print TMP " 1\n";
print TMP " 0.0\n";
print TMP "**END OF DALTON INPUT\n";
close TMP;

open(TMP,">IN_DAL_POLAR_MP2_S");
print TMP "**DALTON INPUT\n";
print TMP ".RUN RESPONSE\n";
print TMP "**WAVE FUNCTIONS\n";
print TMP ".HF\n";
print TMP ".MP2\n";
print TMP "**RESPONSE\n";
print TMP ".SOPPA\n";
print TMP "*LINEAR\n";
print TMP ".DIPLEN\n";
print TMP ".FREQUENCIES\n";
print TMP " 1\n";
print TMP " 0.0\n";
print TMP "**END OF DALTON INPUT\n";
close TMP;
}


sub dodaltondisp {
# make IN_DALTON
open(TMP,">IN_DALTON");
print TMP "basis=6-311G**\n";
print TMP "ATOMBASIS\n";
print TMP "molecule\n";
print TMP "Disp calculation, no charge field\n";
close TMP;
# make IN_DISP
open(TMP,">IN_DISP");
print TMP "**DALTON INPUT\n";
print TMP ".RUN RESPONSE\n";
print TMP "**WAVE FUNCTIONS\n";
print TMP ".HF\n";
print TMP "**RESPONSE\n";
print TMP ".MAXRM\n";
print TMP "400\n";
print TMP "*C6\n";
print TMP ".MAXMOM\n";
print TMP "20\n";
print TMP ".GSLEGN\n";
print TMP ".DIPLEN\n";
print TMP ".PRINT\n";
print TMP "10\n";
print TMP "**END OF DALTON INPUT\n";
close TMP;

$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
$n++;
}
close $fh;
$n=$n-1;

$qmprog=$store[1];
for (my $ic=0; $ic<=$n; $ic++) {
 if ("$store[$ic]" eq "Is long range dispersion accounted for?")  {
  $ie=$ic+1;
  $disp=$store[$ie];
 }
}

if("$disp" eq "Y") {

&SMFAdaldispinputs_MAC;
#&SMFArundispdal;
}

}


sub SMFAdaldispinputs_MAC  {

open(IN,">>INLISTDAL");
 my $nfrgs = `ls grp.*.mol | wc -l`;
  for (my $ic=1; $ic<=$nfrgs; $ic++)
  {
   my $filebase="grp.$ic";
   &Daltonify("$filebase.mol");
   system ("ln -s IN_DISP grp.$ic.dal");
   print IN "$filebase.mol\n";
  }
close IN;
}

sub Daltonify  {

 $filename="$_[0]";

$pwd=`pwd`;
chomp $pwd;

$header="IN_DALTON";

open (my $comfile, "$filename") or die "Unable to open $filename";
@coords=<$comfile>;
close $comfile;

open (my $comfile, ">$filename") or die "Unable to open $filename for writing";
read_dump_headers($header,$comfile);
if (defined $basis)
{
   dump_coords($comfile,$basis,$H_basis);
}
else
{
   print $comfile @coords;
}
   print $comfile "\n";

close $comfile;
}

sub read_dump_headers 
{
        $Headerfile = $_[0];
             $fh = $_[1];

         if ($filename =~ /mol$/){    #dalton molecule files always end in .mol
             $filename =~ s/\.[^.]*$//;
         }

        open(HEADER,"$Headerfile") or die "Unable to open $Headerfile";
        @HeaderLines = <HEADER>;
        while ($HLine = shift (@HeaderLines))
        {
                if ($HLine =~ /molecule/) #find the title line
                {
                        print $fh "$filename\n";
                }
                elsif($HLine =~ /^H(B|b)asis=/) # stick a line with the basis in this file
                {
                   $H_basis=$';
                   chomp $H_basis;
                }
                elsif($HLine =~ /^(B|b)asis=/) # stick a line with the basis in this file
                {
                   $basis=$';
                   chomp $basis;
                }
                else
                {
                        print $fh "$HLine";
                }
        }
        close HEADER;
        if(!defined($H_basis))
        {
           #If we haven't defined a separate basis set for hydrogen atoms, then hydrogen atoms use the same
           #           #basis as all other atoms
           $H_basis=$basis
        }
}

sub dump_coords
{
             $fh = $_[0];
        $basis=$_[1];
        $H_basis=$_[2];

        while ($HLine = shift (@coords))
        {
           #print $HLine;
           if ($HLine =~ /Basis/) #find the Basis line 
              {
                 $tmp_start_of_line=$`;
                 if($tmp_start_of_line =~ /^Charge=1.0/) #Line refers to hydrogen
                 {
                    print $fh "$tmp_start_of_line Basis=$H_basis\n";
                 }
                 else
                 {
                   print $fh "$tmp_start_of_line Basis=$basis\n";
                }
              }
              else
              {
                 #           print $fh;
                      print $fh "$HLine";
              }
        }

}

sub SMFArundispdal  {

system("rm -f *.dal");

 my $nfrgs = `ls grp.*.mol | wc -l`;
  for (my $ic=1; $ic<=$nfrgs; $ic++)
  {
   system ("ln -s IN_DISP grp.$ic.dal");
   system("dalton -N 1 -o grp.$ic.out grp.$ic");
  }
}

