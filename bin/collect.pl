#!/usr/bin/perl

$deriv =`awk 'NR==2,NR==2 {print \$1}' IN_JOBTYPE`;
chomp($deriv);
$deriv=~ s/^\s+|\s+$//g;

$line=`awk /ncpus=/ pbsfile`;
@bline=split("=",$line);
$ncpus=$bline[1];
chomp($ncpus);
$ncpus=~ s/^\s+|\s+$//g;

$nfrgs = `ls COORD* | wc -l`;
open(FR,">FragDerivatives");
print FR "$deriv\n";
print FR "$nfrgs";
close FR;

for ($i=1;$i <= $ncpus;$i++) {
system("cat FragDerivatives$i >> FragDerivatives");
system("rm FragDerivatives$i");
}

if($deriv >= 2) {
 open(DR,">FragDipderivs");
 print DR "$nfrgs";
 close DR;
 for ($i=1;$i <= $ncpus;$i++) {
 system("cat FragDipderivs$i >> FragDipderivs");
 system("rm FragDipderivs$i");
 }
}



system("ls -1 ab.*.log > abfiles");
$ab=`wc -l abfiles`;
@fields = split / /, $ab;
$nfragsab=$fields[0];

if($nfragsab > 0) {

 if($deriv >= 1) {
 system("echo $nfragsab > abforces");
 for ($i=1;$i <= $ncpus;$i++) {
 system("cat abforces$i >> abforces");
 system("rm abforces$i");
 }
}

 if($deriv >= 2) {
  system("echo $nfragsab > abhessians");
  for ($i=1;$i <= $ncpus;$i++) {
  system("cat abhessians$i >> abhessians");
  }
 }

$Totab=0;
for ($i=1;$i <= $ncpus;$i++) {
 $tot=`awk 'NR==1,NR==1 {print \$1}' NearInt_energy$i`;
 chomp($tot);
 $tot=~ s/^\s+|\s+$//g;
 system("rm NearInt_energy$i");
 $Totab=$Totab+$tot;
 }
 system("echo $Totab > NearInt_energy");
}

