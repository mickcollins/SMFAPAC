#!/usr/bin/perl

$CODEDIR=`awk 'NR==1,NR==1 {print \$1}' CODENAME`;
chomp($CODEDIR);
$CODEDIR=~ s/^\s+|\s+$//g;

$EXEDIR=`awk 'NR==1,NR==1 {print \$1}' EXENAME`;
chomp($EXEDIR);
$EXEDIR=~ s/^\s+|\s+$//g;


&dodaltonpolar;

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
 } else {
   if ($nchgs eq 0) {
&SMFAmkdalpolar;
   }
}

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

sub SMFAmkdalpolar {
system("ls -1 nb.*.0-polar.com > polarfiles");
$filename='polarfiles';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";


while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 system("$EXEDIR/convert2dalton < nb.$ic1.0-polar.com > nb.$ic1.0-polar.mol");
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
}











