#!/usr/bin/perl

$llim=$ARGV[0];
$llim=~ s/^\s+|\s+$//g;
$ulim=$ARGV[1];
$ulim=~ s/^\s+|\s+$//g;

$C="C";
$C="$C$llim";

$nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;

&SMFAmkchinputs_nwc;

sub SMFAmkchinputs_nwc {
if ($llim > 0) {
system("ls -1 charge.*.coord > chargefiles$C");
$filename="chargefiles";
$filename="$filename$C";
open($fh,"<$filename");
while ($row = <$fh>){
if ( $. >= $llim && $. <= $ulim )
{
chomp($row);
  @fields=split(/\./,$row);
  $ic=$fields[1];

  open(TMP,">tmpfile$C");
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
  system("tail -$mf charge.$ic.coord >> tmpfile$C");
  open(TMP,">>tmpfile$C");
  print TMP "END\n";
  print TMP "BASIS\n";
  close TMP;
  system("cat NWCbasis >> tmpfile$C");
  open(TMP,">>tmpfile$C");
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
   $n1=~ s/^\s+|\s+$//g;
   system("cat charge.$n1.grp >> tmpfile$C");
   }
  close CHID;
  open(TMP,">>tmpfile$C");
  print TMP "END\n";
    }
  print TMP "ESP\n";
  print TMP "recalculate\n";
  print TMP "END\n";
  print TMP "TASK SCF ENERGY\n";
  print TMP "TASK ESP";
  close TMP;
  system("mv tmpfile$C charge.$ic.nw");
}
}
# end the llim if
}
#end
}



