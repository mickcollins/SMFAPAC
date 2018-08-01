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

$CODEDIR=`awk 'NR==1,NR==1 {print \$1}' CODENAME`;
chomp($CODEDIR);
$CODEDIR=~ s/^\s+|\s+$//g;

$EXEDIR=`awk 'NR==1,NR==1 {print \$1}' EXENAME`;
chomp($EXEDIR);
$EXEDIR=~ s/^\s+|\s+$//g;


&SMFAgaminputs_MAC;

sub SMFAgaminputs_MAC {

    $nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;
# nchgs = 0 for GAMESS
$nchgs = 0;

$deriv=`awk 'NR==2,NR==2 {print \$1}' IN_JOBTYPE`;
chomp $deriv;
$deriv=~ s/^\s+|\s+$//g;

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
open(JUNK,">gamjunk$F");
print JUNK "$prop\n";
close JUNK;
system("sed 's/RHF/UHF/' gamjunk$F > gamjunk2$F");
open(JUNK2,"<gamjunk2$F");
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
  system("tail -$mf GAMbasis > gbfile$F");

# start with the main fragments
if($llim[1] > 0) {
  my $nfrgs = `ls COORD* | wc -l`;
  chomp $nfrgs;
# for (my $ic=1; $ic<=$nfrgs; $ic++)
  for (my $ic=$llim[1]; $ic<=$ulim[1]; $ic++)
  {
  $ch=`awk 'NR==1,NR==1{print \$1}' COORD$ic`;
  $mult=`awk 'NR==1,NR==1{print \$3}' COORD$ic`;
  chomp $ch;
  chomp $mult;
  $mult=~ s/^\s+|\s+$//g;
  open(TMP,">tmpfile$F");
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
    system("grep SYSTEM GAMextras >> tmpfile$F");
   }else{
   system("cat GAMextras >> tmpfile$F");
   }
  }

  open(TMP,">>tmpfile$F");
  print TMP " \$DATA\n";
  print TMP " GAMESS job\n";
  print TMP "c1\n";

  my $lf=`wc -l COORD$ic`;
  chomp $lf;
  my @kf=split(/\s/,"$lf");
  $mf=$kf[0];
  $mf=$mf-1;
  system("tail -$mf COORD$ic > tfile$F");
  system("$EXEDIR/convertcoords4games < tfile$F > gamesscoords$F");
# have to interleave gamesscoords with atombasis and blanl line
  open(TMP1,"<gamesscoords$F");

  while (<TMP1>) {
   my $line1=$_;
   chomp $line1;
   my @bits=split(/\s+/,$line1);
   $ele=$bits[0];
   print TMP "$line1\n";
if ($ans eq "N") {
  open(GAMB,"<gbfile$F");
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
  system("cat gbfile$F >> tmpfile$F");
  open(TMP,">>tmpfile$F");
}
}
  print TMP "\$END\n";
  close TMP;
  close TMP1;
  system("mv tmpfile$F FRAG$ic.inp");
  }
#end the llim[1] if
}
# now the ab.*.* files
if($llim[2] > 0) {
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

  $ch=`awk 'NR==3,NR==3{print \$1}' ab.$ic1.$ic2.com`;
  $mult=`awk 'NR==3,NR==3{print \$2}' ab.$ic1.$ic2.com`;
  chomp $ch;
  chomp $mult;
  $mult=~ s/^\s+|\s+$//g;
  open(TMP,">tmpfile$A");
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
    system("grep SYSTEM GAMextras >> tmpfile$A");
   }else{
   system("cat GAMextras >> tmpfile$A");
   }
  }

  open(TMP,">>tmpfile$A");
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
  system("tail -$mf ab.$ic1.$ic2.com > tfile1$A");
  $mf=$mf-1;
  system("head -$mf tfile1$A > tfile$A");
  system("$EXEDIR/convertcoords4games < tfile$A > gamesscoords$A");
  open(TMP1,"<gamesscoords$A");
  while (<TMP1>) {
   my $line1=$_;
   chomp $line1;
   my @bits=split(/\s+/,$line1);
   $ele=$bits[0];
   print TMP "$line1\n";
if ($ans eq "N") {
  open(GAMB,"<gbfile$F");
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
  system("cat gbfile$F >> tmpfile$A");
  open(TMP,">>tmpfile$A");
}
}
  print TMP "\$END\n";
  close TMP;
  close TMP1;
  system("mv tmpfile$A ab.$ic1.$ic2.inp");
}
  }
close $fh;
# end the llim[2] if
}
# now the nb.* files
if($llim[3] > 0) {
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
  $ch=`awk 'NR==1,NR==1{print \$1}' nb.$ic1.0.com`;
  $mult=`awk 'NR==1,NR==1{print \$2}' nb.$ic1.0.com`;
  chomp $ch;
  chomp $mult;
  $mult=~ s/^\s+|\s+$//g;
  open(TMP,">tmpfile$N");
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
    system("grep SYSTEM GAMextras >> tmpfile$N");
   }else{
   system("cat GAMextras >> tmpfile$N");
   }
  }

  open(TMP,">>tmpfile$N");
  print TMP " \$DATA\n";
  print TMP " GAMESS job\n";
  print TMP "c1\n";


  my $lf=`wc -l nb.$ic1.0.com`;
  chomp $lf;
  my @kf=split(/\s/,"$lf");
  $mf=$kf[0];
  $mf=$mf-1;
  system("tail -$mf nb.$ic1.0.com > tfile1$N");
  $mf=$mf-1;
  system("head -$mf tfile1$N > tfile$N");
  system("$EXEDIR/convertcoords4games < tfile$N > gamesscoords$N");
  open(TMP1,"<gamesscoords$N");
  while (<TMP1>) {
   my $line1=$_;
   chomp $line1;
   my @bits=split(/\s+/,$line1);
   $ele=$bits[0];
   print TMP "$line1\n";
if ($ans eq "N") {
  open(GAMB,"<gbfile$F");
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
  system("cat gbfile$F >> tmpfile$N");
  open(TMP,">>tmpfile$N");
}
}
  print TMP "\$END\n";
  close TMP;
  close TMP1;
  system("echo ' \$stone' >> tmpfile$N");
  system("echo '  bigexp=0.0' >> tmpfile$N");
  system("echo '  atoms' >> tmpfile$N");
  system("echo ' \$end' >> tmpfile$N");
  system("mv tmpfile$N nb.$ic1.0.inp");
}
  }
close $fh;
# end the llim[3] if
}
# the nb polar files are constructed for DALTON
# # as  for NWChem 
# # see rundalpolar

}







