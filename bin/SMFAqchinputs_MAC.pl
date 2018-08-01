#!/usr/bin/perl

for ($i=1;$i <=4;$i++) {
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
$P="P";
$F="$F$llim[1]";
$A="$A$llim[2]";
$N="$N$llim[3]";
$P="$P$llim[4]";

&SMFAqchinputs_MAC;

sub SMFAqchinputs_MAC {


$nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;
$nchgs=~ s/^\s+|\s+$//g;
$deriv=`awk 'NR==2,NR==2 {print \$1}' IN_JOBTYPE`;
chomp $deriv;
$deriv=~ s/^\s+|\s+$//g;

#open(IN,">INLIST");

# start with the main fragments
if ($llim[1] > 0) {
  my $nfrgs = `ls COORD* | wc -l`;
 chomp $nfrgs;

# for (my $ic=1; $ic<=$nfrgs; $ic++)
  for (my $ic=$llim[1]; $ic<=$ulim[1]; $ic++)
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
   $n1=~ s/^\s+|\s+$//g;
   system("cat charge.$n1.grp >> FRAG$ic.com");
   }
   close CHID;
   system("echo '\$end' >> FRAG$ic.com");
   }
# if deriv=2 (hessian) then we need to add a first job
# # for the gradient
if($deriv == 2) {
 system("cp FRAG$ic.com tmpfile1$F");
 system("sed 's/JOBTYPE FREQ/JOBTYPE FORCE/' tmpfile1$F > tmpfile$F");
 system("echo '\@\@\@' >> tmpfile$F");
 system("echo '' >> tmpfile$F");
 system("cat FRAG$ic.com >> tmpfile$F");
 system("mv tmpfile$F FRAG$ic.com");
# print IN "FRAG$ic.com\n";
}

  }
# end the llim[1] if
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
 $ic1=$fields[1];
 $ic2=$fields[2];
 my $ch1=`awk 'NR==1,NR==1 {print \$1}' nb.$ic1.0.com`;
 chomp($ic1);
 my $ch2=`awk 'NR==1,NR==1 {print \$1}' nb.$ic2.0.com`;
 chomp($ic2);
 my $chtot=$ch1 + $ch2;
 my $chm=$ch1*$ch2;
 if($chm != 0 && $chtot == 0) {
  &abcom;
 }else{ 
 system("cp IN_QCH tmpfile$A");
 system("echo '\$molecule' >> tmpfile$A");
# get the sign of the ab nb combination
 $_=`head -1 ab.$ic1.$ic2.com`;
 s/!Isg_coeff=/\$comment!Isg_coeff=/;
 my $comm=$_;
 my $lf=`wc -l ab.$ic1.$ic2.com`;
 chomp $lf;
 my @kf=split(/\s/,"$lf");
 $mf=$kf[0];
 $mf=$mf-1;
 system("head -$mf ab.$ic1.$ic2.com > tfile$A");
 $mf=$mf-2;
 system("tail -$mf tfile$A  >> tmpfile$A");
 system("echo '\$end' >> tmpfile$A");
 system("echo '$comm' >> tmpfile$A");
 }
 if($nchgs >= 1){
  system("echo '\$external_charges' >> tmpfile$A");
  open (CHID,"chab.$ic1.$ic2");
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $n1=~ s/^\s+|\s+$//g;
   system("cat charge.$n1.grp >> tmpfile$A");
   }
  close CHID;
  system("echo '\$end' >> tmpfile$A");
  }
 system("mv tmpfile$A ab.$ic1.$ic2.com");
if($deriv == 2) {
 system("cp ab.$ic1.$ic2.com tmpfile1$A");
 system("sed 's/JOBTYPE FREQ/JOBTYPE FORCE/' tmpfile1$A > tmpfile$A");
 system("echo '\@\@\@' >> tmpfile$A");
 system("echo '' >> tmpfile$A");
 system("cat ab.$ic1.$ic2.com >> tmpfile$A");
 system("mv tmpfile$A ab.$ic1.$ic2.com");
}
# print IN "ab.$ic1.$ic2.com\n";

}
    }
close $fh;
# end the llim[2] if
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
 system("cp IN_QCH tmpfile$N");
 system("echo '\$molecule' >> tmpfile$N");
 my $lf=`wc -l nb.$ic1.0.com`;
 chomp $lf;
 my @kf=split(/\s/,"$lf");
 $mf=$kf[0];
 $mf=$mf-1;
 system("head -$mf nb.$ic1.0.com >> tmpfile$N");
 system("echo '\$end' >> tmpfile$N");
 if($nchgs >= 1){
  system("echo '\$external_charges' >> tmpfile$N");
  open (CHID,"chnb.$ic1");
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $n1=~ s/^\s+|\s+$//g;
   system("cat charge.$n1.grp >> tmpfile$N");
   }
  close CHID;
  system("echo '\$end' >> tmpfile$N");
  }
 system("mv tmpfile$N nb.$ic1.0.com");
if($deriv == 2) {
 system("cp nb.$ic1.0.com tmpfile1$N");
 system("sed 's/JOBTYPE FREQ/JOBTYPE FORCE/' tmpfile1$N > tmpfile$N");
 system("echo '\@\@\@' >> tmpfile$N");
 system("echo '' >> tmpfile$N");
 system("cat nb.$ic1.0.com >> tmpfile$N");
 system("mv tmpfile$N nb.$ic1.0.com");
}
# print IN "nb.$ic1.0.com\n";
}
 }
close $fh;
# end the llim[3] if
}
# now any polarisability jobs
if ($llim[4] > 0) {
 system("ls -1 nb.*.0-polar.com > polarfiles$P");
$filename='polarfiles';
$filename="$filename$P";
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
if ($. >= $llim[4] && $. <= $ulim[4] )
{
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 system("cp IN_QCHPOLAR tmpfile$P");
 system("echo '\$molecule' >> tmpfile$P");
 my $lf=`wc -l nb.$ic1.0-polar.com`;
 chomp $lf;
 my @kf=split(/\s/,"$lf");
 $mf=$kf[0];
 $mf=$mf-1;
 system("head -$mf nb.$ic1.0-polar.com >> tmpfile$P");
 system("echo '\$end' >> tmpfile$P");
 system("mv tmpfile$P nb.$ic1.0-polar.com");
# print IN "nb.$ic1.0-polar.com\n";
}
}
close $fh;
# end the llim[4] if
}
#close IN;
# end of subroutine
}


sub abcom {

my $lf=`wc -l IN_QCH`;
chomp $lf;
 my @kf=split(/\s/,"$lf");
 $mf=$kf[0];
 $mf=$mf-1;
 system("head -$mf IN_QCH > tmpfile$A");
 system("echo 'SCF_GUESS       FRAGMO' >> tmpfile$A");
 system("echo 'SCF_PRINT_FRGM  FALSE >> tmpfile$A");
 system("echo 'MOM_START   1' >> tmpfile$A");
 system("echo '\$end' >> tmpfile$A");
 system("echo '\$molecule' >> tmpfile$A");
# get the sign of the ab nb combination
 $_=`head -1 ab.$ic1.$ic2.com`;
 s/!Isg_coeff=/\$comment!Isg_coeff=/;
 my $comm=$_;
 my $line3=`awk 'NR==3,NR==3 {print \$0}' ab.$ic1.$ic2.com`;
 chomp($line3);
 system("echo $line3 >> tmpfile$A");
 system("echo '--' >> tmpfile$A");
 my $lf=`wc -l nb.$ic1.0.com`;
 chomp($lf);
 my @kf=split(/\s/,"$lf");
 $mf=$kf[0];
 $mf=$mf-1;
 system("head -$mf nb.$ic1.0.com >> tmpfile$A");
 system("echo '--' >> tmpfile$A");
 my $lf=`wc -l nb.$ic2.0.com`;
 chomp($lf);
 my @kf=split(/\s/,"$lf");
 $mf=$kf[0];
 $mf=$mf-1;
 system("head -$mf nb.$ic2.0.com >> tmpfile$A");
 system("echo '\$end' >> tmpfile$A");
 system("echo '$comm' >> tmpfile$A");


}


