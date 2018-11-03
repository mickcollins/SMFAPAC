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
if($disp ne "Y") {
$disp="N";
}
$deriv=0;
if ("$jobtype" eq "SP") {$deriv=0}
if ("$jobtype" eq "FORCE") {$deriv=1}
if ("$jobtype" eq "FREQ") {$deriv=2}
if ("$jobtype" eq "OPT") {$deriv=3}
if ("$jobtype" eq "TS") {$deriv=4}
if ("$jobtype" eq "SCAN") {$deriv=5}
print "\n";
RETRYbasis:
print "QChem allows a different basis set for each element\n";
print "Do you want all elements to have the same basis (Y or N) ? ";
my $ans=<STDIN>;
chomp $ans;
$ans=uc($ans);
if($ans ne "Y" && $ans ne "N") {
 print "The answer must be Y or N\n";
 system("sleep 2");
 goto RETRYbasis;
}
if($ans ne "N") {$ans="Y"};
if ($ans eq "Y") {
print "Enter the basis for all atoms\n";
$basis=<STDIN>;
chomp $basis;
system("rm -f QCHbasis");
}else{
print "Enter the basis (check availability in the QChem library)\n";
print "for the following elements\n";
open (TMP,"<labels0");
system("rm -f QCHbasis");
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
#print TMP "METHOD $method\n";
print TMP "METHOD HF\n";
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

# make a special IN_QCH for step0 in an opt/ts
open TMP,">IN_QCH_step0";
print TMP "\$rem\n";
print TMP "JOBTYPE FREQ\n";
print TMP "IPRINT 10000000\n";
print TMP "METHOD HF\n";
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
system("cat QCHbasis >> IN_QCH_step0");
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
$_=~ s/^\s+|\s+$//g;
   my @secline = split (/\s+/, $_);
   $secen = "$secline[3]";
  }
   if (/Total energy in the final/) {
$_=~ s/^\s+|\s+$//g;
    my @totline = split (/\s+/, $_);
    $ent = "$totline[8]";
    $TotEn = $ent - $secen;
   }
   if(/correlation energy/) {
     if(/Total/) {
$_=~ s/^\s+|\s+$//g;
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

#open(IN,">INLIST");

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
}
# print IN "FRAG$ic.com\n";

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
# print IN "nb.$ic1.0-polar.com\n";
}
close $fh;
#close IN;
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

#$filename='ABSETTINGS';
#open($fh,'<:encoding(UTF-8)',$filename)
#or die "could not open '$filename' $!";
#$n=0;
#while ($row = <$fh>){
#chomp $row;
#$store[$n]=$row;
#$n++;
#}
#close $fh;
#$n=$n-1;
#$nproc0=$store[$n];

$nproc0=1;

#for (my $ic=0; $ic<=$n; $ic++) {
# if ("$store[$ic]" eq "Is long range dispersion accounted for?")  {
#  $ie=$ic+1;
#  $disp=$store[$ie];
# }
#}

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
 print CHK "$stem.fchk\n";
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
#system ("ln -s IN_DISP $stem.dal");
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




