#!/usr/bin/perl

&getscanenergies;


sub getscanenergies {

open(OUT1,">junk");
open(TMP,"<OUT_SMFA");
 my $skip=0;
 my $writit=0;
 $icount=0;
 while (<TMP>) {
  if (/The final geometry has a total energy of/) {
   $skip=$. + 1;
   $icount=$icount+1;
  }
  if ($. == $skip) {$writit=1};
  if( $writit == 1) {
   $energy=$_;
   chomp $energy;
   print OUT1 "$icount      $energy\n";
   $writit=0;
  }
 }
close OUT1;
close TMP;
open(TMP,">>OUT_SMFA");
print TMP "\n";
print TMP "Summary of Scan Energies\n";
print TMP "\n";
close TMP;
system("cat junk >> OUT_SMFA");
}


