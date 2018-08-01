#!/usr/bin/perl

$CODEDIR=`awk 'NR==1,NR==1 {print \$1}' CODENAME`;
chomp($CODEDIR);
$CODEDIR=~ s/^\s+|\s+$//g;

$EXEDIR=`awk 'NR==1,NR==1 {print \$1}' EXENAME`;
chomp($EXEDIR);
$EXEDIR=~ s/^\s+|\s+$//g;

&extractch_nwc;

sub extractch_nwc  {

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

system("tail -$Nat charge.$ic.esp > tmpfile");
      open (NPA_IN,"<tmpfile");
      open (NPA_OUT,">charge.$ic.npa");
      print NPA_OUT "charge.$ic.coord\n";
      print NPA_OUT "$Nat\n";
      $n=0;
      while(<NPA_IN>) {
      $line="$_";
      @bits=split(/\s+/,"$line");
      $q=$bits[4];
      print NPA_OUT "$x[$n] "," $y[$n] "," $z[$n] "," $q\n";
      $n++;
}
      close NPA_OUT;
      system("$EXEDIR/adjustcaponly < charge.$ic.npa");
}

#end
}

