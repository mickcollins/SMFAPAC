#!/usr/bin/perl

for ($i=1;$i <=3;$i++) {
$i1=2*($i-1);
$i2=$i1+1;
$llim[$i]=$ARGV[$i1];
$llim[$i]=~ s/^\s+|\s+$//g;
$ulim[$i]=$ARGV[$i2];
$ulim[$i]=~ s/^\s+|\s+$//g;
}

$F="F";
$A="A";
$N="N";
$F="$F$llim[1]";
$A="$A$llim[2]";
$N="$N$llim[3]";

&SMFAnwcinputs_MAC;


sub SMFAnwcinputs_MAC {


$nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;
$nchgs=~ s/^\s+|\s+$//g;
$deriv=`awk 'NR==2,NR==2 {print \$1}' IN_JOBTYPE`;
chomp $deriv;
$deriv=~ s/^\s+|\s+$//g;

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

$deriv =`awk 'NR==2,NR==2 {print \$1}' IN_JOBTYPE`;
chomp $deriv;

# start with the main fragments
if ($llim[1] > 0) {
  my $nfrgs = `ls COORD* | wc -l`;
  chomp $nfrgs;
# for (my $ic=1; $ic<=$nfrgs; $ic++)
  for (my $ic=$llim[1]; $ic<=$ulim[1]; $ic++)
  {
  open(TMP,">tmpfile$F");
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
  system("tail -$mf COORD$ic >> tmpfile$F");
  open(TMP,">>tmpfile$F");
  print TMP "END\n";
if($nchgs >= 1){
  open (CHID,"chFRAG.$ic");
  print TMP "bq\n";
  close TMP;
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $n1=~ s/^\s+|\s+$//g;
   system("cat charge.$n1.grp >> tmpfile$F");
   }
  close CHID;
  open(TMP,">>tmpfile$F");
    }
  if($nchgs >= 1){print TMP "END\n"};

  print TMP "BASIS\n";
  close TMP;
  system("cat NWCbasis >> tmpfile$F");
  open(TMP,">>tmpfile$F");
  print TMP "END\n";
 open(FH,"<NWCcommands");
 while (<FH>) {
  if( $mult == 2 ) {
   s/SINGLET/$mform/;
  }
  print TMP "$_";
 }
  close TMP;
  system("mv tmpfile$F FRAG$ic.nw");
}
# end the $llim[1] if
}
# now the ab.*.* files
if ($llim[2] > 0) {
 system("ls -1 ab.*.com > abfiles$A");
$filename='abfiles';
$filename="$filename$A";
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
if ($. >= $llim[2] && $. <= $ulim[2] )
{
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 my $ic2=$fields[2];
  open(TMP,">tmpfile$A");
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
 system("head -$mf ab.$ic1.$ic2.com > tfile$A");
 $mf=$mf-3;
 system("tail -$mf tfile$A  >> tmpfile$A");
 open(TMP,">>tmpfile$A");
 print TMP "END\n";
if($nchgs >= 1){
  open (CHID,"chab.$ic1.$ic2");
  print TMP "bq\n";
  close TMP;
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $n1=~ s/^\s+|\s+$//g;
   system("cat charge.$n1.grp >> tmpfile$A");
   }
  close CHID;
  open(TMP,">>tmpfile$A");
    }
  if($nchgs >= 1){print TMP "END\n"};
  print TMP "BASIS\n";
  close TMP;
  system("cat NWCbasis >> tmpfile$A");
  open(TMP,">>tmpfile$A");
  print TMP "END\n";
 open(FH,"<NWCcommands");
 while (<FH>) {
  if( $mult == 2 ) {
   s/SINGLET/$mform/;
  }
  print TMP "$_";
 }
  close TMP;
  system("mv tmpfile$A ab.$ic1.$ic2.nw");
}
 }
close $fh;
# end the $llim[2] if
}
# now the nb.* files
if ($llim[3] > 0) {
system("ls -1 nb.*.0.com > nbfiles$N");
my $pwd=`pwd`;
chomp $pwd;
$filename='nbfiles';
$filename="$filename$N";
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
if ($. >= $llim[3] && $. <= $ulim[3] )
{
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
  open(TMP,">tmpfile$N");
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
 system("head -$mf nb.$ic1.0.com > tfile$N");
 $mf=$mf-1;
 system("tail -$mf tfile$N >> tmpfile$N");
 open(TMP,">>tmpfile$N");
 print TMP "END\n";
if($nchgs >= 1){
  open (CHID,"chnb.$ic1");
  print TMP "bq\n";
  close TMP;
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $n1=~ s/^\s+|\s+$//g;
   system("cat charge.$n1.grp >> tmpfile$N");
   }
  close CHID;
  open(TMP,">>tmpfile$N");
    }
  if($nchgs >= 1){print TMP "END\n"};
  print TMP "BASIS\n";
  close TMP;
  system("cat NWCbasis >> tmpfile$N");
  open(TMP,">>tmpfile$N");
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
  system("mv tmpfile$N nb.$ic1.0.nw");
}
}
close $fh;
# end the llim[3] if
}
# the nb polar files are constructed for DALTON
# # as NWChem does not have polarizabilities
# # see rundalpolar
}






