#!/usr/bin/perl

$llim=$ARGV[0];
$llim=~ s/^\s+|\s+$//g;
$ulim=$ARGV[1];
$ulim=~ s/^\s+|\s+$//g;
$C="C";
$C="$C$llim";

&SMFAmkchinputs_gau;



sub SMFAmkchinputs_gau {

$nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;
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
  }
# end the llim if
}
# end sub
}

 sub ltrim { my $s = shift; $s =~ s/^\s+//;       return $s };

