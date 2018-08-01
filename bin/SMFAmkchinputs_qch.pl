#!/usr/bin/perl

$llim=$ARGV[0];
$llim=~ s/^\s+|\s+$//g;
$ulim=$ARGV[1];
$ulim=~ s/^\s+|\s+$//g;
$C="C";
$C="$C$llim";

&SMFAmkchinputs_qch;

sub SMFAmkchinputs_qch {
if($llim > 0) {
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
  system("cp IN_QCHNPA charge.$ic.com");
  system("echo '\$molecule' >> charge.$ic.com");
  system("cat charge.$ic.coord >> charge.$ic.com");
  system("echo '\$end' >> charge.$ic.com");
  if ( -s "chL0.$ic")  {
   system("echo '\$external_charges' >> charge.$ic.com");
   open (CHID,"chL0.$ic");
   foreach my $n1 (<CHID>) {
   chomp $n1;
   $n1=~ s/^\s+|\s+$//g;
   system("cat charge.$n1.grp >> charge.$ic.com");
   }
   close CHID;
   system("echo '\$end' >> charge.$ic.com");
  }
}
}
# end the llim if
}
# end sub
}





