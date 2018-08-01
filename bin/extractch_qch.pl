#!/usr/bin/perl

$CODEDIR=`awk 'NR==1,NR==1 {print \$1}' CODENAME`;
chomp($CODEDIR);
$CODEDIR=~ s/^\s+|\s+$//g;

$EXEDIR=`awk 'NR==1,NR==1 {print \$1}' EXENAME`;
chomp($EXEDIR);
$EXEDIR=~ s/^\s+|\s+$//g;


&extractch_qch;

sub extractch_qch {
      my @jobs = `ls charge.*.coord`;
      foreach my $chg_coord (@jobs) {
      @fields=split(/\./,$chg_coord);
      $ic=$fields[1];

      open (CHGFIL,"$chg_coord");
      $n=0;
      while (<CHGFIL>) {
      chomp;
      if ($. == 1) {
      ($dum,$chg,$mlt) = split /\s+/, $_;
      } else {
      ($symb[$n],$x[$n],$y[$n],$z[$n])=split(/\s+/,$_);
      $n++;
      }
      $Nat = $n;
      }
      close CHGFIL;


$start=0;
$end=0;

      open (NPA_IN,"charge.$ic.log");
      $n=0;
      while (<NPA_IN>) {
      if (/Summary of Natural Population Analysis/){
      $start=$. + 5;
      $end=$start+$Nat;
      }
      if (($. > $start) && ($. <= $end)) {
      ($dum,$dum,$dum,$chgarray[$n]) = split (/\s+/, $_);
      $n++;
      }
      }
      close NPA_IN;
      open (NPA_OUT,">charge.$ic.npa");
      print NPA_OUT "charge.$ic.coord\n";
      print NPA_OUT "$Nat\n";

      for ($n=0; $n < $Nat; $n++) {
      print NPA_OUT "$x[$n] $y[$n] $z[$n] $chgarray[$n]\n";
      }
      close NPA_OUT;

      system("$EXEDIR/adjustcaponly < charge.$ic.npa");

      }
# end sub
}




