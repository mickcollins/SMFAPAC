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

&SMFAgauinputs_MAC;

sub SMFAgauinputs_MAC  {

$nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;
$nchgs=~ s/^\s+|\s+$//g;
$deriv=`awk 'NR==2,NR==2 {print \$1}' IN_JOBTYPE`;
chomp $deriv;
$deriv=~ s/^\s+|\s+$//g;

my $pwd=`pwd`;
chomp $pwd;

#open(IN,">INLIST");
# start with the main fragments
if ($llim[1] > 0) {
  my $nfrgs = `ls COORD* | wc -l`;
 chomp $nfrgs;
# for (my $ic=1; $ic<=$nfrgs; $ic++)
  for (my $ic=$llim[1]; $ic<=$ulim[1]; $ic++)
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
  system("cp IN_G09 IN_G09_temp$F");
  if($ffile[0] ne "0") {
   system("sed 's/nosymm/charge nosymm/' IN_G09_temp$F > IN_G09_current$F");
  }else{
   system("cp IN_G09 IN_G09_current$F");
  }
  if($deriv >= 2) {
  system("echo '%chk=$pwd/FRAG$ic.chk' > FRAG$ic.com");
  }
  system("cat IN_G09_current$F >> FRAG$ic.com");
  system("echo 'Bonded fragment', $ic, ' with any necessary charges\n' >> FRAG$ic.com");
  system("cat $my_coord >> FRAG$ic.com");
  system("echo ' ' >> FRAG$ic.com");
    if($nchgs >= 1){
  open (CHID,"chFRAG.$ic");
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $m1 = &ltrim($n1);
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
#print IN "FRAG$ic.com\n";
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
# find out if there are any external charges
  $ffile[0]=0;
  if( -s "chab.$ic1.$ic2") {
  $efile=`wc -l chab.$ic1.$ic2`;
  chomp $efile;
  @ffile=split(/\s+/,$efile);
  }
# if necessary add charge to the G09 header
  system("cp IN_G09 IN_G09_temp$A");
  if($ffile[0] ne "0") {
   system("sed 's/nosymm/charge nosymm/' IN_G09_temp$A > IN_G09_current$A");
  }else{
   system("cp IN_G09 IN_G09_current$A");
  }

 system("cp IN_G09_current$A tmpfile$A");
 system("cat ab.$ic1.$ic2.com >> tmpfile$A");
    if($nchgs >= 1){
  open (CHID,"chab.$ic1.$ic2");
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $m1 = &ltrim($n1);
   system("cat charge.$m1.grp >> tmpfile$A");
   }
  close CHID;
  system("echo ' ' >> tmpfile$A");
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
open (TMP,">>tmpfile$A");
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
 system("mv tmpfile$A ab.$ic1.$ic2.com");
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
 system("echo '%chk=$pwd/nb.$ic1.0.chk' > tmpfile$N");
# find out if there are any external charges
  $ffile[0]=0;
  if( -s "chnb.$ic1") {
  $efile=`wc -l chnb.$ic1`;
  chomp $efile;
  @ffile=split(/\s+/,$efile);
  }
# if necessary add charge to the G09 header
  system("cp IN_G09 IN_G09_temp$N");
  if($ffile[0] ne "0") {
   system("sed 's/nosymm/charge nosymm/' IN_G09_temp$N > IN_G09_current$N");
  }else{
   system("cp IN_G09 IN_G09_current$N");
  }
  system("cat IN_G09_current$N >> tmpfile$N");

 system("echo 'nonbonded monomer with any necessary charges\n' >> tmpfile$N");
 system("cat nb.$ic1.0.com >> tmpfile$N");
    if($nchgs >= 1){
  open (CHID,"chnb.$ic1");
  foreach my $n1 (<CHID>) {
   chomp $n1;
   $m1 = &ltrim($n1);
   system("cat charge.$m1.grp >> tmpfile$N");
   }
  close CHID;
  system("echo ' ' >> tmpfile$N");
    }
# now add the gen basis if necessary
  if ($ans ne 0) {
# get only the basis for the elements present
open (TMP,">>tmpfile$N");
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
 system("mv tmpfile$N nb.$ic1.0.com");
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
 system("cp IN_POLAR tmpfile$P");
 system("echo 'polarisability job for nonbonded monomer\n' >> tmpfile$P");
 system("cat nb.$ic1.0-polar.com >> tmpfile$P");
# now add the gen basis if necessary
  if ($ans ne 0) {
open (TMP,">>tmpfile$P");
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
 system("mv tmpfile$P nb.$ic1.0-polar.com");
# print IN "nb.$ic1.0-polar.com\n";
}
}
close $fh;
# end the $llim[4] if
}
#close IN;
# end of subroutine
}

sub ltrim { my $s = shift; $s =~ s/^\s+//;       return $s };





