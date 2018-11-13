#!/usr/bin/perl

$CODEDIR=`awk 'NR==1,NR==1 {print \$1}' CODENAME`;
chomp($CODEDIR);
$CODEDIR=~ s/^\s+|\s+$//g;

$EXEDIR=`awk 'NR==1,NR==1 {print \$1}' EXENAME`;
chomp($EXEDIR);
$EXEDIR=~ s/^\s+|\s+$//g;

&fragonly;

sub fragonly {

system("rm charge.*");
system("rm OUT_NCOORD_FRAG_Lev0");
system("rm OUT_FRAGMENT*");
system("rm OUT_L1L1_data");
system("rm CHARGECOORDS");
system("rm OUT_CHARGEDGROUPS");
system("rm VIEWCHARGEDGROUPS");
system("rm OUT_NCOORD_FRAG");
system("rm fragnu*");
system("rm frags.ou*");
system("rm signs.ou*");
system("rm *COORD*");
system("rm chL*");
system("rm chFR*");
system("rm chnb*");
system("rm chab*");
system("rm grp.*");
system("rm nb.*");
system("rm ab*");
system("rm ABbg*");
system("rm FRAG*");
system("rm Frag*");
system("rm Lev0*");
system("rm Lev1*");
system("rm El*");
system("rm Near*");
system("rm Ind*");
system("rm Disp*");
system("rm Doub*");
system("rm -f IN_HBONDS");
system("rm -f list.*");
system("rm -f tmpf*");
system("rm -f tfil*");
system("rm -f IN_G09_temp*");
system("rm -f IN_G09_current*");
system("rm -f polarfiles*");
system("rm -f OUT_ABGROUPS");
system("rm -f STEPNUMBER");
system("rm -f gaminput*");
system("rm -f gauinput*");
system("rm -f gauchinput*");
system("rm -f nwcinput*");
system("rm -f nwcchinput*");
system("rm -f qchinput*");
system("rm -f qchchinput*");
system("rm -f runchg*");
system("rm -f runpar*");
system("rm -f rundal*");
system("rm -f TIMELOADING");

  system("echo 0 > ONEONLY");

  open(TMP2,">IN_FRACT");
  print TMP2 "Enter the Level of fragmentation\n";
  print TMP2 "0\n";
  print TMP2 "Enter the cutoff value for nonbonded interactions\n";
  print TMP2 "1.1\n";
  close TMP2;
  system("$EXEDIR/cutsmfa > rub0");
  system("mv OUT_NCOORD_FRAG OUT_NCOORD_FRAG_Lev0");
  system("mv OUT_FRAGMENTS OUT_FRAGMENTS_Lev0");
  system("mv fragnum fragnum_Lev0");
  system("mv frags.out frags.out_Lev0");
  system("mv signs.out signs.out_Lev0");
  my @coords = `ls COORD*`;
   foreach my $my_coord (@coords) {
    chomp $my_coord;
    system("mv $my_coord Lev0_${my_coord}");
   }

  system("$EXEDIR/findHandCbonds_VdW");

  system("echo 1 > ONEONLY");

  open(TMP2,">IN_FRACT");
  print TMP2 "Enter the Level of fragmentation\n";
  print TMP2 "1\n";
  print TMP2 "Enter the cutoff value for nonbonded interactions\n";
  print TMP2 "1.1\n";
  close TMP2;
  system("$EXEDIR/cutsmfa > rub1");
  system("mv OUT_NCOORD_FRAG OUT_NCOORD_FRAG_Lev1");
  system("mv OUT_FRAGMENTS OUT_FRAGMENTS_Lev1");
  system("mv fragnum fragnum_Lev1");
  system("mv frags.out frags.out_Lev1");
  system("mv signs.out signs.out_Lev1");
  system("rm -f chFRAG*");
  my @coords1 = `ls COORD*`;
   foreach my $my_coord (@coords1) {
    chomp $my_coord;
    system("mv $my_coord Lev1_${my_coord}");
   }

  system("echo 2 > ONEONLY");

  system("cp IN_FRACT_LevX IN_FRACT");
  system("$EXEDIR/cutsmfa > rubX");
  system("mv OUT_FRAGMENTS OUT_FRAGMENTS_LevX");

$package =`awk 'NR==1,NR==1 {print \$0}' IN_PACKAGE`;
chomp $package;
if($package == 1) {
open(CHG,">OUT_CHARGEDGROUPS");
print CHG " For GAMESS the number of groups with formal charges is zero\n";
print CHG "0\n";
close CHG;
}

&run_L1L1MAC;

&catELECTS;

# don't re-run setprocs
system("rm -f INLIST*");
system("$EXEDIR/orderlist");

if ( -s "INLISTCHG" ) {
 system("cp INLISTCHG INLISTCHG_original");
}
if ( -s "INLIST" ) {
 system("cp INLIST INLIST_original");
}
if ( -s "INLISTDAL" ) {
 system("cp INLISTDAL INLISTDAL_original");
}


open(TMP2,">READY2GO");
print TMP2 "2\n";
close(TMP2);

      return 1;
}

    sub run_L1L1MAC { #ok

     system("rm -f fort.23");
     system("rm -f fort.24");
     system("rm -f fort.25");
     system("rm -f fort.26");
     system("$EXEDIR/L1L1_SMFA > rubL1L1");
     if( -s "fort.23") {
     system("cat fort.23 >> OUT_ELECTRONS");
     }
     if( -s "fort.24") {
     system("cat fort.24 >> OUT_ELECTRONS");
     }
     if( -s "fort.25") {
     system("cat fort.25 >> OUT_ELECTRONS");
     }
     if( -s "fort.26") {
     system("cat fort.26 >> OUT_ELECTRONS");
     }


     system("stty sane");
     return 1;
    }

sub catELECTS {

$already=`grep 'Total numbers of atoms and electrons' OUT_SMFA  | wc -l`;
chomp $already;
if ( $already == 0 ) {
system("cat OUT_ELECTRONS_SUMMARY >> OUT_SMFA");
#&largest;
}else{
 $Celects=`awk 'NR==1,NR==1 {print \$8}' OUT_ELECTRONS_SUMMARY`;
 chomp $Celects;
 $Cfrags=`awk 'NR==3,NR==3 {print \$1}' OUT_ELECTRONS_SUMMARY`;
 chomp $Cfrags;
 open(TMP,"<OUT_SMFA");
 while(<TMP>) {
  if(/Total numbers of atoms and electrons/) {
   $line=$_;
  }
 }
 close TMP;
 @num=split(/\s+/,$line);
 $length= scalar(@num);
 $Oelects=$num[$length-1];
 if ($Celects != $Oelects) {
 system("cat OUT_ELECTRONS_SUMMARY >> OUT_SMFA");
 open(TMP,">>OUT_SMFA");
 print TMP "The number of electrons in the molecule has changed !!!!\n";
 print TMP "This has occurred because SMFA has introduced new (spurious) charges in the molecule\n";
 print TMP "The new charges are listed in the file IN_CHARGES,\n";
 print TMP "created by the program Preparegeom\n";
 print TMP "Such spurious charges may have come about because Preparegeom detected new bonds or\n";
 print TMP "broken bonds compared to the original geometry (during an optimization or scan)\n";
 print TMP "Check the current geometry in the file Newcoords.xyz\n";
 print TMP "Perhaps the chosen scan path is unreasonable?\n";
 print TMP "If the current geometry seems OK, you can avoid this error by specifying\n";
 print TMP "the correct charge of the atoms concerned in the input, and re-running the\n";
 print TMP "optimisation/scan.\n";
 print TMP "SMFA is aborted by this error.\n";
 close TMP;
 exit(0);
 }
 open(TMP,"<OUT_SMFA");
 $skip=0;
 while(<TMP>) {
 if (/The number of fragments is/) {
  $skip=$. + 1;
 }
 if ($. == $skip) {
  $line=$_;
 }
 }
 close TMP;
 chomp $line;
 $line=~ s/^\s+|\s+$//g;
 $Ofrags=$line;
 print "$Ofrags   $Cfrags\n";
 if ( $Ofrags != $Cfrags) {
  open(TMP,">>OUT_SMFA");
  print TMP "The fragmentation has changed from the original set of fragments.\n";
  print TMP "Presumeably, the bonding has changed during an optimization or scan.\n";
  print TMP "This is not a fatal error. However, it means that the energy will not change\n";
  print TMP "continuously during the optimization/scan.\n";
  print TMP "If you want to avoid such a discontinuity, then\n";
  print TMP "1. Compare the original geometry with the current geometry in Newcoords.xyz\n";
  print TMP "2. Set the 'unusual bonding' in the input to ensure the bonding is the same\n";
  print TMP "   throughout the optimisation/scan, and\n";
  print TMP "3. Repeat the optimization/scan.\n";
 }
}

}









