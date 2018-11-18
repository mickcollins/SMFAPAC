######################### view the input #######################A

sub viewinp  { 
system("stty sane");
system("clear");
system("stty echo");
print "Waiting\n";
open(SMP1,">OUT_SMFA");
print SMP1 "              INPUT to SMFA\n";
print SMP1 "\n";
open(TMP,"<IN_JOBTYPE");
while (<TMP>) {
if ($. == 2) {$deriv="$_"};
}
close TMP;
system("cp IN_JOBTYPE IN_JOBTYPE_original");
chomp $deriv;
if ($deriv == 0) {
print SMP1 "The purpose of the calculation is to evaluate the energy\n";
print SMP1 "\n";
}
if ($deriv == 1) {
print SMP1 "The purpose of the calculation is to evaluate the energy gradient\n";
print SMP1 "\n";
}
if ($deriv == 2) {
print SMP1 "The purpose of the calculation is to evaluate the frequencies\n";
print SMP1 "\n";
}
if ($deriv == 3) {
print SMP1 "The purpose of the calculation is to find a structure with minimum energy\n";
print SMP1 "\n";
}
if ($deriv == 4) {
print SMP1 "The purpose of the calculation is to find a saddle point geometry (TS)\n";
print SMP1 "\n";
}
if ($deriv == 5) {
print SMP1 "The purpose of the calculation is to perform an optimised scan\n";
print SMP1 "\n";
}
if ($deriv == 3) {
close SMP1;
system("cat OPTSTEPS >> OUT_SMFA");
$junk=`cat IN_CONSTRAINTS >> OUT_SMFA`;
open(SMP1,">>OUT_SMFA");
print SMP1 "\n";
}
if ($deriv == 4) {
close SMP1;
system("cat OPTSTEPS >> OUT_SMFA");
$junk=`cat IN_CONSTRAINTS >> OUT_SMFA`;
system("cat TSATOMS >> OUT_SMFA");
open(SMP1,">>OUT_SMFA");
print SMP1 "\n";
}
if ($deriv == 5) {
close SMP1;
system("cat IN_OPTSCAN >> OUT_SMFA");
system("cat OPTSTEPS >> OUT_SMFA");
open(SMP1,">>OUT_SMFA");
print SMP1 "\n";
}
$filename='IN_FRACT_LevX';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
chomp $row;
print SMP1 "$row\n";
}
close $fh;
print SMP1 "\n";
$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
if ($n eq 0) {print SMP1 "$store[0]\n";}
if ($n == 1) {
if ($store[1] == 4) {print SMP1 "QChem\n";}
if ($store[1] == 3) {print SMP1 "NWChem\n";}
if ($store[1] == 2) {print SMP1 "GAUSSIAN09\n";}
if ($store[1] == 1) {print SMP1 "GAMESS US\n";}
}
if ($n > 1) {print SMP1 "$row\n";}
$n++;
}
print SMP1 "\n";
$filename='xyzFILENAME';
print SMP1 "The molecular coordinates are in file\n";
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while ($row = <$fh>){
chomp $row;
if($n == 0) {print SMP1 "$row\n";}
if($n == 1) {
print SMP1 "The number of atoms in the molecule is\n";
print SMP1 "$row\n";
}
$n++;
}
close $fh;
print SMP1 "\n";
$filename='IN_PRESERVE_BONDS';
if (-s $filename) {
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
chomp $row;
print SMP1 "$row\n";
}
close $fh;
print SMP1 "\n";
}
$filename='IN_MINUSBONDS';
if (-s $filename) {
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
chomp $row;
print SMP1 "$row\n";
}
close $fh;
print SMP1 "\n";
}
$filename='IN_DOUBLE';
if (-s $filename) {
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
chomp $row;
print SMP1 "$row\n";
}
close $fh;
print SMP1 "\n";
}
$filename='IN_SINGLE';
if (-s $filename) {
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
chomp $row;
print SMP1 "$row\n";
}
close $fh;
print SMP1 "\n";
}
$filename='IN_EXTRABONDS';
if (-s $filename) {
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
chomp $row;
print SMP1 "$row\n";
}
close $fh;
print SMP1 "\n";
}
$filename='IN_SPECIFIED_CHARGES';
if (-s $filename) {
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
chomp $row;
print SMP1 "$row\n";
}
close $fh;
print SMP1 "\n";
}



system "rm -f OUT_METALGROUPS";
system "rm -f OUT_NONMETALCHARGES";
system("rm -f WARNINGS");

&callPrepare;

if ( -s "IN_SOLCHARGES" ) {
system("$EXEDIR/SOLCH");
print SMP1 "Solvent molecules will contribute embedded charges to all calculations,\n";
print SMP1 "as well as appearing explicitly in ab initio calculations\n";
print SMP1 "\n";
close SMP1;
system("cat IN_SOLCHARGES >> OUT_SMFA");
}else{
close SMP1;
}

open(SMP1,">>OUT_SMFA");

print SMP1 "\n";
print SMP1 "             EXAMINATION of the MOLECULAR COORDINATES shows:\n";
print SMP1 "\n";

$filename='molcharge';
if (-s $filename) {
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
chomp $row;
print SMP1 "$row\n";
}
close $fh;
print SMP1 "\n";
}

close SMP1;
system("cat -s WARNINGS >> OUT_SMFA");
open(SMP1,">>OUT_SMFA");
print SMP1 "\n";

$filename='OUT_METALGROUPS';
if (-s $filename) {
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
chomp $row;
print SMP1 "$row\n";
}
close $fh;
print SMP1 "\n";
}
$filename='OUT_NONMETALCHARGES';
if (-s $filename) {
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
chomp $row;
print SMP1 "$row\n";
}
close $fh;
print SMP1 "\n";
}

print "Done\n";
print "View and check the input in file OUT_SMFA\n";
print "Hit RETURN to continue ";

# note input processed in READY2GO;
open(TMP2,">READY2GO");
print TMP2 "1\n";
close(TMP2);

my $rub=<STDIN>;
      return 1;

}

######################### get cartesian coordinates #############################

    sub get_coords { #ok
       &endwin();
      system("stty sane");
      SMFAgeominp();
      return 1;
    }

######################### get the ab method data #############################

    sub get_abdata { #ok
       &endwin();
      system("stty sane");
      SMFAabinp();
      return 1;
    }

sub SMFAabinp {
system("clear");
system("stty echo");

print "What quantum chemistry program package will you use?\n";
print "The current choices are:\n";
print "1 GAMESS US\n";
print "2 Gaussian09\n";
print "3 NWChem\n";
print "4 QChem\n";
print "Choose a number\n";
$package=<STDIN>;
if ($package eq "\n") { return 1 };
chomp $package;
if($package != 1 && $package != 2 && $package != 3 && $package != 4) {
print "INCORRECT PROGRAM PACKAGE NUMBER, SMFA WILL ABORT\n";
system("sleep 5");
}
      open TMP, ">IN_PACKAGE";
      print TMP "$package\n";
      close TMP;
if ( $package == 1 ) {&readgamparams};
if ( $package == 2 ) {&readgauparams};
if ( $package == 3 ) {&readnwcparams};
if ( $package == 4 ) {&readqchparams};

open(TMP2,">READY2GO");
print TMP2 "0\n";
close(TMP2);
}

######################### get bonding information #############################

    sub get_bonds { #ok
      &endwin();
      system("stty sane");
      SMFAbondinp();
      return 1;
     }

######################### get charge information #############################

    sub get_charges { #ok
       &endwin();
      system("stty sane");
      SMFAchargeinp();
      return 1;
     }

######################### get fragmentation info  #############################

     sub get_frag { #ok
            &endwin();
           system("stty sane");
           SMFAfraginp();
           return 1;
     }

######################### Fragment molecule only ########################
  sub fragonly {
system("stty sane");
system("clear");
system("stty echo");
# check if input has been processed;
$filename='READY2GO';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($ready = <$fh>){
chomp $ready;
$ready2=trim($ready);
}
close $fh;

if ($ready2 < "1") {&viewinp};

system("rm -f charge.*");
system("rm -f OUT_NCOORD_FRAG_Lev0");
system("rm -f OUT_FRAGMENT*");
system("rm -f OUT_L1L1_data");
system("rm -f CHARGECOORDS");
system("rm -f OUT_CHARGEDGROUPS");
system("rm -f VIEWCHARGEDGROUPS");
system("rm -f OUT_NCOORD_FRAG");
system("rm -f fragnu*");
system("rm -f frags.ou*");
system("rm -f signs.ou*");
system("rm -f *COORD*");
system("rm -f chL*");
system("rm -f chFR*");
system("rm -f chnb*");
system("rm -f chab*");
system("rm -f grp.*");
system("rm -f nb*");
system("rm -f ab*");
system("rm -f ABbg*");
system("rm -f FRAG*");
system("rm -f Frag*");
system("rm -f Lev0*");
system("rm -f Lev1*");
system("rm -f El*");
system("rm -f Near*");
system("rm -f Ind*");
system("rm -f Disp*");
system("rm -f Doub*");
system("rm -f INLIST*");
system("rm -f list.*");
system("rm -f tmpf*");
system("rm -f tfil*");
system("rm -f IN_G09_temp*");
system("rm -f IN_G09_current*");
system("rm -f polarfiles*");
system("rm -f OUT_ABGROUPS");
system("rm -f ONEONLY");
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

print "Waiting for fragmentation at Level = 0\n";
  system("rm -f IN_HBONDS");

  open(TMP2,">IN_FRACT");
  print TMP2 "Enter the Level of fragmentation\n";
  print TMP2 "0\n";
  print TMP2 "Enter the cutoff value for nonbonded interactions\n";
  print TMP2 "1.1\n";
  close TMP2;

  system("echo 0 > ONEONLY");

  system("$EXEDIR/cutsmfa > rub0");
  $error=`grep "Fragmentation will freeze" rub0 | wc -l`;
  if ($error > 0) {
  print "SMFA cannot fragment this molecule\n";
  print "SMFA will abort\n";
  system("echo 'SMFA cannot fragment this molecule' >> OUT_SMFA");
  system("echo 'SMFA will abort' >> OUT_SMFA");
  exit(0);
  }
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

print "Waiting for fragmentation at Level = 1\n";

  system("$EXEDIR/findHandCbonds_VdW");
  system("cat -s IN_HBONDS >> OUT_SMFA");

  open(TMP2,">IN_FRACT");
  print TMP2 "Enter the Level of fragmentation\n";
  print TMP2 "1\n";
  print TMP2 "Enter the cutoff value for nonbonded interactions\n";
  print TMP2 "1.1\n";
  close TMP2;

  system("echo 1 > ONEONLY");

  system("$EXEDIR/cutsmfa > rub1");
  $error=`grep "Fragmentation will freeze" rub1 | wc -l`;
  if ($error > 0) {
  print "SMFA cannot fragment this molecule\n";
  print "SMFA will abort\n";
  system("echo 'SMFA cannot fragment this molecule' >> OUT_SMFA");
  system("echo 'SMFA will abort' >> OUT_SMFA");
  exit(0);
  }

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

print "Waiting for fragmentation at the requested Level\n";
  system("cp IN_FRACT_LevX IN_FRACT");

  system("echo 2 > ONEONLY");

  system("$EXEDIR/cutsmfa > rubX");
  $error=`grep "Fragmentation will freeze" rubX | wc -l`;
  if ($error > 0) {
  print "SMFA cannot fragment this molecule\n";
  print "SMFA will abort\n";
  system("echo 'SMFA cannot fragment this molecule' >> OUT_SMFA");
  system("echo 'SMFA will abort' >> OUT_SMFA");
  exit(0);
  }

  system("mv OUT_FRAGMENTS OUT_FRAGMENTS_LevX");
  &viewfrags;
# get the number of fragments
$nfragsX=`awk 'NR==2,NR==2 {print \$0}' frags.out`;
chomp($nfragsX);
$nfragsX=~ s/^\s+|\s+$//g;

open(SMP2,">>OUT_SMFA");
print SMP2 "\n";
print SMP2 "             FRAGMENTATION RESULTS\n";
print SMP2 "\n";
print SMP2 "A description of the fragments can be found in file frags.out\n";
print SMP2 "The coordinates of the fragments in xyz format are contained\n";
print SMP2 "in seefrags (in suitable format for graphics).\n";
print SMP2 "\n";

if( $nfragsX == 1 ) {
print SMP2 "                   WARNING\n";
print SMP2 "\n";
print SMP2 "THIS MOLECULE DOES NOT FRAGMENT AT THIS VALUE OF LEVEL\n";
print SMP2 "                SMFA is ABORTING\n";
print "\n";
print "                   WARNING\n";
print "\n";
print "THIS MOLECULE DOES NOT FRAGMENT AT THIS VALUE OF LEVEL\n";
print "                SMFA is ABORTING\n";
system("sleep 5");
exit(0);
}

close SMP2;

# GAMESS does not use charged groups as embedded charges, so before L1L1_smfa
# we have to put the number of charges to zero
$package =`awk 'NR==1,NR==1 {print \$0}' IN_PACKAGE`;
chomp $package;
if($package == 1) {
open(CHG,">OUT_CHARGEDGROUPS");
print CHG " For GAMESS the number of groups with formal charges is zero\n";
print CHG "0\n";
close CHG;
}

print "Waiting for creation of the non-bonded ab intio calculations\n";
&run_L1L1MAC;

&catELECTS;

system("$EXEDIR/setprocs >> OUT_SMFA");
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

print "Done\n";
print "See OUT_SMFA for information about the fragments.\n";
print "\n";
print "OUT_SMFA also contains a recommendation regarding the number of cpus\n";
print "to be used for a parallel execution of the electronic structure calculations.\n";
print "\n";
print "Hit RETURN to continue ";
$rub=<STDIN>;
# record that all the fragmentation steps have been done;
open(TMP2,">READY2GO");
print TMP2 "2\n";
close(TMP2);

      return 1;
}

sub SMFAgeominp {
system("clear");
system("stty echo");

print "Name of the file containing the cartesian coordinates in xyz format - see Users Guide\n";
$file=<STDIN>;
if ($file eq "\n") { return 1 };
chomp $file;

if (!-s $file) {
   print "$file does not exist or is empty.\n";
   system("sleep 5");
   die;
}
# remove leading blank spaces from each line of $file
&trimxyz("$file");

# new geometry needs to remove old input files
system("rm -f IN_DOUBLE IN_EXTRABONDS IN_HBONDS IN_ISOTOPES");
system("rm -f IN_MINUSBONDS IN_PRESERVE_BONDS IN_SINGLE ");
system("rm -f IN_SPECIFIED_CHARGES IN_CHARGES IN_CAPS");
system("rm -f IN_CONSTRAINTS*");


  my ($num_at, $charge, $mass);
  $n = 0;
  open(FL,$file) or die "Cannot open $file!";
  while(<FL>){
     if( $_ eq "\n") {next};
     chomp;
     next if $. < 3;
      if (/END/) {next};
     $n++;
     ($symb[$n],$x[$n],$y[$n],$z[$n])=split(/\s+/,$_);
  }
  $Nat = $n;

open(TMP1,">molxyz");
print TMP1 "$Nat\n";
print TMP1 "$file\n";
for ($n=1; $n <= $Nat; $n++) {
     print TMP1 "$symb[$n],$x[$n],$y[$n],$z[$n]\n";
}
close(TMP1);
open(TMP1,">xyzFILENAME");
print TMP1 "$file\n";
print TMP1 "$Nat\n";
close(TMP1);

open(TMP2,">READY2GO");
print TMP2 "0\n";
close(TMP2);

system("stty sane");
}

sub SMFAbondinp {

system("clear");
system("stty echo");
print "If you want to accept all SMFA defaults regarding bond definitions,\n";
print "simply hit the RETURN key to skip this input\n";
print "\n";
RETRYdouble:
print "Do you want to specify some bonds as multiple bonds,\n";
print "that would normally be taken as single bonds (Y/N or hit RET to skip all)?\n";
  $ifdouble=<STDIN>;
if ($ifdouble eq "\n") {return 1}
chomp $ifdouble;
$ifdouble=uc($ifdouble);
if($ifdouble ne "Y" && $ifdouble ne "N") {
 print "The answer must be Y or N\n";
 system("sleep 2");
 goto RETRYdouble;
}
$numberd=0;
if ($ifdouble eq "Y") {
print "Enter the number of (unusual) user-specified double bonds (usually 0)\n";
  $numberd=<STDIN>;
chomp $numberd;
print "Enter the atom numbers for these double bonded atoms (eg 27   243)\n";
for ($i=1; $i <= $numberd; $i++) {
$line=<STDIN>;
chomp $line;
($ndb[$i],$mdb[$i])= split /\s+/, $line;
    }
  open(TMP2,">IN_DOUBLE");
  print TMP2 "The number of specified multiple bonds is\n";
  print TMP2 "$numberd\n";
  print TMP2 "The atoms connected by these bonds are\n";
  for ($i=1; $i <= $numberd; $i++) {
       print TMP2 "$ndb[$i] $mdb[$i]\n";
    }
  close(TMP2);
} else {
  open(TMP2,">IN_DOUBLE");
  print TMP2 "The number of specified multiple bonds is\n";
  print TMP2 "$numberd\n";
  close(TMP2);
}

RETRYsingle:
print "Do you want to specify some bonds as single bonds,\n";
print "that would normally be taken as multiple bonds (Y/N)?\n";
  $ifsingle=<STDIN>;
chomp $ifsingle;
$ifsingle=uc($ifsingle);
if ($ifsingle ne "Y" && $ifsingle ne "N") {
 print "The answer must be Y or N\n";
 system("sleep 2");
 goto RETRYsingle;
}
$numbers=0;
if($ifsingle eq "Y") {
print "Enter the number of (unusual) user-specified single bonds (usually 0)\n";
  $numbers=<STDIN>;
chomp $numbers;
print "Enter the atom numbers for these single bonded atoms (eg 86   97)\n";
for ($i=1; $i <= $numbers; $i++) {
$line=<STDIN>;
chomp $line;
($ndb[$i],$mdb[$i])= split /\s+/, $line;
           }
  open(TMP2,">IN_SINGLE");
  print TMP2 "The number of specified single bonds is\n";
  print TMP2 "$numbers\n";
  print TMP2 "The atoms connected by these bonds are\n";
  for ($i=1; $i <= $numbers; $i++) {
       print TMP2 "$ndb[$i] $mdb[$i]\n";
    }
  close(TMP2);
} else {
  open(TMP2,">IN_SINGLE");
  print TMP2 "The number of specified single bonds is\n";
  print TMP2 "$numbers\n";
  close(TMP2);
}

RETRYextra:
print "Do you want to specify that single bonds exist\n";
print "that would normally not be considered as bonds (Y/N)?\n";
  $ifbonds=<STDIN>;
chomp $ifbonds;
$ifbonds=uc($ifbonds);
if($ifbonds ne "Y" && $ifbonds ne "N") {
 print "The answer must be Y or N\n";
 system("sleep 2");
 goto RETRYextra;
}
$numberd=0;
if ($ifbonds eq "Y"){
print "Enter the number of (unusual) user-specfied extra bonds (usually 0)\n";
  $numberd=<STDIN>;
chomp $numberd;
print "Enter the atom numbers for the atoms connected by these extra bonds (eg 147   15)\n";
for ($i=1; $i <= $numberd; $i++) {
$line=<STDIN>;
chomp $line;
($ndb[$i],$mdb[$i])= split /\s+/, $line;
    }
  open(TMP3,">IN_EXTRABONDS");
  print TMP3 "The number of extra single bonds is\n";
  print TMP3 "$numberd\n";
  print TMP3 "The atom numbers for the connected atoms are\n";
  for ($i=1; $i <= $numberd; $i++) {
       print TMP3 "$ndb[$i] $mdb[$i]\n";
    }
  close(TMP3);
} else {
  open(TMP3,">IN_EXTRABONDS");
  print TMP3 "The number of extra single bonds is\n";
  print TMP3 "$numberd\n";
  close(TMP3);
}
RETRYamide:
print "By default SMFA assumes that the amide moiety is treated as a single group\n";
print "Do you want to treat the Amide CN bond as a single bond only (Y/N)?\n";
  $ifamide=<STDIN>;
chomp $ifamide;
$ifamide=uc($ifamide);
if($ifamide ne "Y" && $ifamide ne "N") {
 print "The answer must be Y or N\n";
 system("sleep 2");
 goto RETRYamide;
}
$lamide="T";
if ($ifamide eq "Y") {$lamide="F"};
if ($ifamide eq "N") {$lamide="T"};

RETRYcx:
print "By default SMFA assumes that CX2 or CX3 (X=F,Cl,..) is a single group.\n";
print "Do you want to treat these CX bonds as single bonds only (Y/N)?\n";
  $ifCF=<STDIN>;
chomp $ifCF;
$ifCF=uc($ifCF);
if($ifCF ne "Y" && $ifCF ne "N") {
 print "The answer must be Y or N\n";
 system("sleep 2");
 goto RETRYcx;
}
$lcf="T";
if ($ifCF eq "Y") {$lcf="F"};
if ($ifCF eq "N") {$lcf="T"};

  open(TMP4,">IN_PRESERVE_BONDS");
  print TMP4 "Classes of moieties that are treated as united groups\n";
  print TMP4 "$lamide   amides\n";
  print TMP4 "$lcf   CX2 and CX3 (X=F,Cl,..)\n";
  close(TMP4);

  $hbonds="true";
  $ibonds="true";
  $factor="1.0";
RETRYhb:
print "By default SMFA assumes that hydrogen bonds will be treated as regular bonds\n";
print "Do you want to ignore hydrogen bonds (Y/N)?\n";
  $ifhb=<STDIN>;
chomp $ifhb;
$ifhb=uc($ifhb);
if($ifhb ne "Y" && $ifhb ne "N") {
 print "The answer must be Y or N\n";
 system("sleep 2");
 goto RETRYhb;
}
if ($ifhb eq "Y") {$hbonds="false"};
RETRYch:
print "By default SMFA assumes that if two charged groups are close together\n";
print "then they should be considered to be bonded\n";
print "Do you want to ignore this default (Y/N)?\n";
  $ifch=<STDIN>;
chomp $ifch;
$ifch=uc($ifch);
$ifch=~ s/^\s+|\s+$//g;
if($ifch ne "Y" && $ifch ne "N") {
 print "The answer must be Y or N\n";
 system("sleep 2");
 goto RETRYch;
}
if ($ifch eq "Y") {$ibonds="false"};
if ($ifch eq "N") {
  print "Two charged groups are considered close together if the distance\n";
  print "between them is less than the sum of the Van der Waals radii of\n";
  print "the closest atoms multiplied by a factor\n";
  print "Enter the value of this factor (the recommended value is 1.0)\n";
  $factor=<STDIN>;
  chomp $factor;
}
  open(TMP5,">IN_MINUSBONDS");
  print TMP5 "Hydrogen bonds are treated as actual bonds\n";
  print TMP5 "$hbonds\n";
  print TMP5 "Close charged groups are taken as bonded\n";
  print TMP5 "$ibonds\n";
  print TMP5 "The factor which multiplies the sum of the vdW radii\n";
  print TMP5 "$factor\n";
  close(TMP5);
open(TMP2,">READY2GO");
print TMP2 "0\n";
close(TMP2);
system("stty sane");
}

sub SMFAchargeinp {

system("clear");
system("stty echo");
RETRYmetals:
print "Are there metals or other atoms whose charge must be specified\n";
print "or modified (Y/N, or hit RETURN to skip)?\n";
  $ifch=<STDIN>;
if ($ifch eq "\n") { return 1 }
chomp $ifch;
$ifch=uc($ifch);
if($ifch ne "Y" && $ifch ne "N") {
 print "The answer must be Y or N or hit RETURN to skip\n";
 system("sleep 2");
 goto RETRYmetals;
}
if ($ifch eq "Y") {
print "Enter the number of metal/other atoms\n";
  $numberc=<STDIN>;
chomp $numberc;
print "Enter the atom number(s) and associated charge(s), eq 34  2 (ie the 34th atom has a charge of 2)\n";
for ($i=1; $i <= $numberc; $i++) {
$line=<STDIN>;
chomp $line;
($an[$i],$ac[$i])= split /\s+/, $line;
    }
 open(TMP2,">IN_SPECIFIED_CHARGES");
 print TMP2 "The number of metal or unusual atoms is \n $numberc \n";
 print TMP2 "The atom number(s) and associated charge(s)\n";
 for ($i=1; $i <= $numberc; $i++) {
      print TMP2 "$an[$i] $ac[$i]\n";
 }
 close(TMP2);

}
if ($ifch eq "N") {
system ("rm -f IN_SPECIFIED_CHARGES");
}
open(TMP2,">READY2GO");
print TMP2 "0\n";
close(TMP2);
system("stty sane");
}

sub SMFAfraginp {
system("clear");
system("stty echo");

print "Enter the Level of fragmentation\n";
$fraglvl=<STDIN>;
if ($fraglvl eq "\n") {return 1};
chomp $fraglvl;
print "The fragmentation level is $fraglvl \n";
print "Enter the cutoff value for nonbonded interactions\n";
print "For fragmentation level = 2, the recommended value is 0.\n";
print "For fragmentation level > 2, the recommended value is 1.1 - see Users Guide\n";
$cutoff=<STDIN>;
chomp $cutoff;
print "The cutoff level is $cutoff \n";

open(TMP2,">IN_FRACT_LevX");
print TMP2 "The Level of fragmentation is\n";
print TMP2 "$fraglvl\n";
print TMP2 "The cutoff value for nonbonded interactions is\n";
print TMP2 "$cutoff\n";
close(TMP2);

open(TMP2,">READY2GO");
print TMP2 "0\n";
close(TMP2);

system("stty sane");
}

sub runall {

# make pbsfile_bare which has no instructions for loading quant chem packages
# and make LOADQCP (loads base package) and LOADDP (loads dalton)

&makeqmdal;
system("chmod +x LOADQCP");
system("chmod +x LOADDP");

    while (1) {  # always true
     &menu_init($numbered_flag, "SMFA Control ",0,
                ,"\n-submenu: Run Options");
     &menu_item ("Run all calculations sequentially","runallseq");
     &menu_item ("Run all calculations in parallel","runallpar1n");
     &menu_item ("Run all calculations in parallel with multiple nodes","runallparmany");
     &menu_item ("Exit run options","Exit");

     $sel = &menu_display("",$menu_default_row,$menu_default_top,
                             $menu_default_col);

     if ($sel eq "%EMPTY%" ) { die "not enough screen lines \n";}
     if ($sel eq "Exit") {last;}

       if ($sel ne "%UP%" )
        { # call the corresponding subprogram (assuming existence)
         $ret=&$sel();
           if ($ret)
            {
            $window=&initscr();
            &menu_curses_application($window);
            } # end of: if ($ret)
        } # end of: if ($sel ne "%UP%" )
    } # while (1)
  } # end of: sub inp

sub runallpar1n {
system("echo 'N' > MULTINODE");
&runallpar;
}

sub runallparmany {
system("echo 'Y' > MULTINODE");
if ( -s "MULTINODEDATA" ) {
}else{
system("clear");
system("stty sane");
system("cp $CODEDIR/MULTINODEDATA MULTINODEDATA");
print "Before you run a multinode calculation in this directory for the first time,\n";
print "you must edit the MULTINODEDATA file in this directory, as directed by the\n";
print "notes in that file. Alternatively, you can copy a previously edited file from\n";
print "another directory into this directory.\n";
print "\n";
print "Hit RETURN to exit\n";
$ans=<STDIN>;
exit(0);
}
&runallpar;

}

sub runallseq {
system("clear");
system("stty sane");

# check if fragmentation has been done
 $filename='READY2GO';
 open($fh,'<:encoding(UTF-8)',$filename)
 or die "could not open '$filename' $!";
 while ($ready = <$fh>){
 chomp $ready;
 $ready2=&trim($ready);
 }
 close $fh;

 if ($ready2 < "2") {
 print "Before a job can be submitted (or re-submitted), SMFA must check the input\n";
 print "and perform the fragmentation.\n";
 print "SMFA will perform those preliminary steps now.\n";
 system("sleep 5");
 &fragonly;
}
 system("echo 0 > READY2GO");

$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
$n++;
}
$nlast=$n-1;
close $fh;
open(TMP,"<IN_JOBTYPE");
while (<TMP>) {
if ($. == 2) {$deriv="$_"};
}
close TMP;
chomp $deriv;

$qmprog=$store[1];
$nproc0=$store[$nlast];

# abandon sequential optimisations and scans for GAMESS and NWChem
# in favour of parallel with 1 cpu, as at 280418
if( $qmprog == 1 || $qmprog == 3 ) {
 if( $deriv >= 3 ) {
 &runallpar;
 exit(0);
 }
}

if( ($deriv == 3) || ($deriv == 4)) {
# this is an optimisation job
# make the script for submission
system("rm -f OUTLIST*");
system("rm optcode.pl");
system("cp $CODEDIR/SMFA_optscript.pl optcode.pl");
system("cat $CODEDIR/SMFA_subs.pl >> optcode.pl");
system("cat $CODEDIR/SMFA_gam.pl >> optcode.pl");
system("cat $CODEDIR/SMFA_gau.pl >> optcode.pl");
system("cat $CODEDIR/SMFA_nwc.pl >> optcode.pl");
system("cat $CODEDIR/SMFA_qch.pl >> optcode.pl");
system("cat $CODEDIR/SMFA_dal.pl >> optcode.pl");
system("chmod +x optcode.pl");
if($qmprog == 1) {
system("cp pbsfile subrunoptgam");
system("echo './optcode.pl' >> subrunoptgam");
system("chmod +x subrunoptgam");
system("qsub subrunoptgam");
}
if($qmprog == 2) {
system("cp pbsfile subrunoptgau");
system("echo './optcode.pl' >> subrunoptgau");
system("chmod +x subrunoptgau");
system("qsub subrunoptgau");
}
if($qmprog == 3) {
system("cp pbsfile subrunoptnwc");
system("echo './optcode.pl' >> subrunoptnwc");
system("chmod +x subrunoptnwc");
system("qsub subrunoptnwc");
}
if($qmprog == 4) {
system("cp pbsfile subrunoptqch");
system("echo './optcode.pl' >> subrunoptqch");
system("chmod +x subrunoptqch");
system("qsub subrunoptqch");
}
print "\n";
print "Ab initio optimisation job submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}

if($deriv == 5) {
# this is a scan job
# make the script for submission
system("rm -f OUTLIST*");
system("rm optcode.pl");
system("cp $CODEDIR/SMFA_scanscript.pl optcode.pl");
system("cat $CODEDIR/SMFA_subs.pl >> optcode.pl");
system("cat $CODEDIR/SMFA_gam.pl >> optcode.pl");
system("cat $CODEDIR/SMFA_gau.pl >> optcode.pl");
system("cat $CODEDIR/SMFA_nwc.pl >> optcode.pl");
system("cat $CODEDIR/SMFA_qch.pl >> optcode.pl");
system("cat $CODEDIR/SMFA_dal.pl >> optcode.pl");
system("chmod +x optcode.pl");
if($qmprog == 1) {
system("cp pbsfile subrunoptgam");
system("echo './optcode.pl' >> subrunoptgam");
system("chmod +x subrunoptgam");
system("qsub subrunoptgam");
}
if($qmprog == 2) {
system("cp pbsfile subrunoptgau");
system("echo './optcode.pl' >> subrunoptgau");
system("chmod +x subrunoptgau");
system("qsub subrunoptgau");
}
if($qmprog == 3) {
system("cp pbsfile subrunoptnwc");
system("echo './optcode.pl' >> subrunoptnwc");
system("chmod +x subrunoptnwc");
system("qsub subrunoptnwc");
}
if($qmprog == 4) {
system("cp pbsfile subrunoptqch");
system("echo './optcode.pl' >> subrunoptqch");
system("chmod +x subrunoptqch");
system("qsub subrunoptqch");
}
print "\n";
print "Ab initio scan job submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}


if ($qmprog == 1) {
if ($cont eq "RESTART") {
system("cp $CODEDIR/restartgamheader runallgam");
}else{
system("cp $CODEDIR/rungamheader runallgam");
system("rm -f OUTLIST*");
}
system("cat $CODEDIR/SMFA_gam.pl >> runallgam");
system("cat $CODEDIR/SMFA_subs.pl >> runallgam");

if ($cont eq "RESTART") {
system("cp $CODEDIR/restartdalheader runalldal");
}else{
system("cp $CODEDIR/rundalheader runalldal");
}
system("cat $CODEDIR/SMFA_gam.pl >> runalldal");
system("cat $CODEDIR/SMFA_dal.pl >> runalldal");
system("cat $CODEDIR/SMFA_subs.pl >> runalldal");

system("chmod +x runallgam");
system("chmod +x runalldal");
#system("cp pbsfile subrunallgam");
system("cp pbsfile_bare subrunallgam");
system("cat LOADQCP >> subrunallgam");
system("echo './runallgam' >> subrunallgam");

#system("echo 'module unload gamess/2016-08-R1' >> subrunallgam");
#system("echo 'module unload openmpi/1.10.2' >> subrunallgam");
#system("echo 'module load dalton' >> subrunallgam");
system("cat LOADDP >> subrunallgam");
system("echo './runalldal' >> subrunallgam");
system("chmod +x subrunallgam");
system("qsub subrunallgam");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}
if ($qmprog == 2) {
if ($cont eq "RESTART") {
system("cp $CODEDIR/restartgauheader runallgau");
}else{
system("cp $CODEDIR/rungauheader runallgau");
system("rm -f OUTLIST*");
}
system("cat $CODEDIR/SMFA_gau.pl >> runallgau");
system ("cat $CODEDIR/SMFA_dal.pl >> runallgau");
system("cat $CODEDIR/SMFA_subs.pl >> runallgau");
system("chmod +x runallgau");
system("cp pbsfile subrunallgau");
system("echo './runallgau' >> subrunallgau");
system("chmod +x subrunallgau");
system("qsub subrunallgau");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}
if ($qmprog == 3) {
if ($cont eq "RESTART") {
system("cp $CODEDIR/restartnwcheader runallnwc");
}else{
system("cp $CODEDIR/runnwcheader runallnwc");
system("rm -f OUTLIST*");
}
system("cat $CODEDIR/SMFA_nwc.pl >> runallnwc");
system("cat $CODEDIR/SMFA_subs.pl >> runallnwc");

system("cp $CODEDIR/rundalheader runalldal");
system("cat $CODEDIR/SMFA_nwc.pl >> runalldal");
system("cat $CODEDIR/SMFA_dal.pl >> runalldal");
system("cat $CODEDIR/SMFA_subs.pl >> runalldal");

system("chmod +x runallnwc");
system("chmod +x runalldal");
#system("cp pbsfile subrunallnwc");
system("cp pbsfile_bare subrunallnwc");
system("cat LOADQCP >> subrunallnwc");
system("echo './runallnwc' >> subrunallnwc");
#system("echo 'module unload nwchem/6.6' >> subrunallnwc");
#system("echo 'module load dalton' >> subrunallnwc");
system("cat LOADDP >> subrunallnwc");
system("echo './runalldal' >> subrunallnwc");
system("chmod +x subrunallnwc");
system("qsub subrunallnwc");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}
if ($qmprog == 4) {
if ($cont eq "RESTART") {
system("cp $CODEDIR/restartqchheader runallqch");
}else{
system("cp $CODEDIR/runqchheader runallqch");
system("rm -f OUTLIST*");
}
system("cat $CODEDIR/SMFA_qch.pl >> runallqch");
system("cat $CODEDIR/SMFA_dal.pl >> runallqch");
system("cat $CODEDIR/SMFA_subs.pl >> runallqch");
system("chmod +x runallqch");
system("cp pbsfile subrunallqch");
system("echo './runallqch' >> subrunallqch");
system("chmod +x subrunallqch");
system("qsub subrunallqch");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}

open(TMP2,">READY2GO");
print TMP2 "1\n";
close(TMP2);

return 1;
}

# this sub is superceded
sub Lev0_chg_MAC_iter {

system("clear");

$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
print "$store[$n]\n";
$n++;
}

$qmprog=$store[1];


$iter=3;
$count=1;
while ($count <= $iter) {

if ($qmprog == 1) {&SMFAdalnpa_MAC};
if ($qmprog == 2) {&SMFAgamnpa_MAC};
if ($qmprog == 3) {&SMFAgaunpa_MAC};
if ($qmprog == 4) {&SMFAqchnpa_Mac};
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

if ($qmprog == 4) {$ext="log"};
if ($qmprog == 3) {$ext="log"};
if ($qmprog == 2) {$ext="out"};
if ($qmprog == 1) {$ext="log"};

$start=0;
$end=0;

      open (NPA_IN,"charge.$ic.$ext");
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
      $count++;
}

      return 1;
    }

# this sub is superceded
sub SMFAabinitioinputs_MAC {

$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
$n++;
}

$qmprog=$store[1];

if ($qmprog == 3) {&SMFAgauinputs_MAC};
if ($qmprog == 4) {&SMFAqchinputs_MAC};
#if ($qmprog == 2) {dalinp};
####if ($qmprog == 3) {gamnpa};
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


sub anal_MAC {

system("clear");
system("stty sane");

      $nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
      $deriv =`awk 'NR==2,NR==2 {print \$1}' IN_JOBTYPE`;
      $fraglvl =`awk 'NR==2,NR==2 {print \$1}' IN_FRACT`;

chomp $nchgs;
chomp $deriv;
chomp $fraglvl;

system("rm -f combined* ");

$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
$n++;
}
$n=$n-1;
close $fh;

$qmprog=$store[1];
for (my $ic=0; $ic<=$n; $ic++) {
 if ("$store[$ic]" eq "Is long range dispersion accounted for?")  {
  $ie=$ic+1;
  $disp=$store[$ie];
 }
}

if ($qmprog == 1) {

&extract_gam;

&makenbgdma_gam;

# always read polar for GAMESS
 &ReadPolar_dal;
}


if ($qmprog == 2) {

 &extract_gau;

 &makenbgdma_gau;

if ("$disp" eq "Y") {
 &ReadPolar_gau;
 } else {
   if ($nchgs eq 0) {
    &ReadPolar_gau;
   }
}
}

if ($qmprog == 3) {

&extract_nwc;

&makenbgdma_nwc;

if ("$disp" eq "Y") {
 &ReadPolar_dal;
 } else {
   if ($nchgs eq 0) {
    &ReadPolar_dal;
   }
}
}

if ($qmprog == 4) {

&extract_qch;

&makenbgdma_qch;

if ("$disp" eq "Y") {
 &ReadPolar_qch;
 } else {
   if ($nchgs eq 0) {
    &ReadPolar_qch;
   }
}
}

if ( ("$nchgs" eq "0") or ( $qmprog == 1 ) ) {
 system("$EXEDIR/Induct_Mac");
}

 if($disp eq "Y") {
  &extract_disp_dal;
  system("$EXEDIR/Disp");
 }



system("$EXEDIR/Derivs_smfa < FragDerivatives > rubd");
if($deriv >= 1) {
system("$EXEDIR/ABNBderivatives");
}
system("$EXEDIR/Elect_Mac > rube");
if ( ($nchgs > 1) && ( $qmprog != 1)) {
system("$EXEDIR/DoubleC_Mac > rubD");
}
system("$EXEDIR/combineallderivatives");
system("cat SMFA.out.temp >> OUT_SMFA");
system("rm SMFA.out.temp");

if($deriv == 1) {
system("echo '' >> OUT_SMFA");
system("echo ' The energy gradients are contained in the file combinedderivs' >> OUT_SMFA");
}

if($deriv == 2) {
system("$EXEDIR/DipoleDerivs < FragDipderivs > combDipDerivs");
system("rm -f IN_ISOTOPES");
system("$EXEDIR/Freq");
system("echo '' >> OUT_SMFA");
system("echo 'The normal mode frequencies are printed in file FREQUENCIES'  >> OUT_SMFA");
system("echo '' >> OUT_SMFA");
system("echo 'The normal mode eigenvectors are printed in file NORMAL_MODES' >> OUT_SMFA");
system("echo '' >> OUT_SMFA");
system("echo 'A simulated IR spectrum is printed in file SPECTRUM'  >> OUT_SMFA");
system("echo '' >> OUT_SMFA");
system("grep 'The zero-point energy' FREQUENCIES >> OUT_SMFA");
}
system("stty sane");
#system(sleep 5);
}

# superceded sub maybe not
sub make_gdma
{
        my $fh=$_[0];
        print $fh "Title \"$prefix\"\n";

#Check if we need a density keyword
        foreach $GLine (@GDMALines)
        {
           if ($GLine =~ /Density/) #Density line needs to go before File if it's there
           {
              print $fh "$GLine";
           }
        }
        print $fh "File $prefix.fchk\n\n";
        print $fh "Angstrom\n\n";
        foreach $GLine (@GDMALines)
        {
                if ($GLine =~ /Punch/) #find the checkpoint file line
                {
                        chomp $GLine;
                        print $fh " Punch $prefix.punch\n";
                }
                elsif ($GLine =~ /ADD/)#origin line if there is one
                {
                        print $fh @Names;
                }
                elsif ($GLine =~ /Density/) #Density line needs to go before File if it's there
                {
                   #Do nothing

                }
                else
                {
                        print $fh "$GLine";
                }
        }
}


sub solventcharges {

system("stty sane");
system("clear");
RETRYembed:
print "Very polar solvent molecules can induce polarisation in both\n";
print "solute and other solvent molecules.\n";
print "SMFA can evaluate the energy of polar solvents and solutes, at a\n";
print "lower Level of Fragmentation, if embedded charges are used to describe\n";
print "the solvent environment. In order to identify the solvent, SMFA needs the\n";
print "chemical composition of the solvent.\n";
print "If your system contains polar solvent molecules,\n";
print "do you want to specify the use of embedded charges? (Y/N)\n";
$anssol=<STDIN>;
if($anssol eq "\n") {
 system("rm IN_SOLCHARGES");
 return;
}
chomp $anssol;
$anssol=uc($anssol);
if($anssol ne "Y" && $anssol ne "N") {
 print "The answer must be Y or N\n";
 system("sleep 2");
 goto RETRYembed;
}
if($anssol ne "Y") {
 system("rm -f IN_SOLCHARGES");
 return;
}
print "Enter the number of atoms in the solvent molecule\n";
$numsol=<STDIN>;
chomp $numsol;
print "You have to enter each element in the solvent. For example, for water\n";
print "you would enter\n";
print "H \n";
print "H \n";
print "O \n";
print "Now enter the elemental symbols for the atoms in the solvent:\n";
for (my $ij=0;$ij <= $numsol;$ij++ )  {
$ele=<STDIN>;
chomp $ele;
$solatom[$ij]=$ele;
}
print "Waiting\n";
open(TMP,">IN_SOLCHARGES");
print TMP "The number of atoms in the solvent molecule is\n";
print TMP "$numsol\n";
print TMP "The elements in the solvent molecule are\n";
for (my $ik=0;$ik <= $numsol;$ik++ )  {
print TMP "$solatom[$ik]\n";
}
close TMP;

}

sub optinput {

# only need the maximum number of iterations
print "The geometry will be optimised in a number of discrete steps (geometry changes).\n";
print "It is prudent to limit the number of steps employed to limit the total cpu time\n";
print "consumed. The optimisation can always be restarted from the last step reached.\n";
print "Enter the maximum number of steps allowed\n";
$maxsteps=<STDIN>;
chomp $maxsteps;
open (STEP,">OPTSTEPS");
print STEP "The maximum number of optimisation steps is\n";
print STEP "$maxsteps\n";
close STEP;
print "\n";
print "The geometry optimisation can be carried out with constraints on bond lengths,\n";
print "valence bond angles, and dihedral angles,if requested.\n";
print "\n";
print "Do you want constraints?\n";
$yescon=<STDIN>;
chomp $yescon;
$yescon=uc($yescon);
if($yescon ne "Y") {
 system("rm -f IN_CONSTRAINTS");
 return;
}
open(CON,">IN_CONSTRAINTS");
print "Enter the number of bond lengths to be constrained\n";
print CON "Enter the number of bond lengths to be constrained\n";
$nbonds=<STDIN>;
chomp $nbonds;
print CON "$nbonds\n";
if ($nbonds >= 1) {
 print "Enter the atom numbers and bond lengths for each constraint\n";
 print CON "Enter the atom numbers and bond lengths for each constraint\n";
 for ($il=1; $il<=$nbonds; $il++) {
 $line=<STDIN>;
 chomp $line;
 print CON "$line\n";
 }
}
print "Enter the number of angles to be constrained\n";
print CON "Enter the number of angles to be constrained\n";
$nangle=<STDIN>;
chomp $nangle;
print CON "$nangle\n";
if ($nangle >= 1) {
 print "Enter the 3 atom numbers and angle in degrees for each angle\n";
 print CON "Enter the 3 atom numbers and angle in degrees for each angle\n";
 for ($il=1; $il<=$nangle; $il++) {
 $line=<STDIN>;
 chomp $line;
 print CON "$line\n";
 }
}
print "Enter the number of dihedral angles to be constrained\n";
print CON "Enter the number of dihedral angles to be constrained\n";
$ndiheds=<STDIN>;
chomp $ndiheds;
print CON "$ndiheds\n";
if ($ndiheds >= 1) {
 print "Enter the 4 atom numbers and dihedral angle in degrees for each angle\n";
 print CON "Enter the 4 atom numbers and dihedral angle in degrees for each angle\n";
 for ($il=1; $il<=$ndiheds; $il++) {
 $line=<STDIN>;
 chomp $line;
 print CON "$line\n";
 }
}
close CON;

}

sub TSinput {
# get the atoms involved in the TS motion
print "\n";
print "Many saddle points may exist on the potential energy surface\n";
print "of large molecules. To help identify the saddle point (TS) of interest,\n";
print "you must enter the atom numbers of a few atoms (say 3), that you expect to\n";
print "change position as the molecule passes through the TS.\n";
print "\n";
print "Enter the number of atoms you will denote as relevant to the TS\n";
$nTSatoms=<STDIN>;
chomp $nTSatoms;
open(TS,">TSATOMS");
print TS "Enter the number of atoms you will denote as relevant to the TS\n";
print TS "$nTSatoms\n";
if($nTSatoms >=1) {
print "Enter each atom number, one per line\n";
print TS "Enter each atom number, one per line\n";
 for ($il=1; $il<=$nTSatoms; $il++) {
 $line=<STDIN>;
 chomp $line;
 print TS "$line\n";
 }
}
close TS;
}

sub scaninput {
print "The geometry will be optimised for a set of constraints.\n";
print "You must specify a set of configurations defined by a path, a sequence of constraints.\n";
print "These constraints apply to bond lengths, angles, and dihedral angles, which are kept fixed\n";
print "while all other degrees of freedom are optimised.\n";
print "\n";
print "See the Manual for details.\n";
print "\n";
open(CON,">IN_OPTSCAN");
print CON "Enter the number of geometries in the scan\n";
print "Enter the number of geometries in the scan\n";
$ngeom=<STDIN>;
chomp $ngeom;
print CON "$ngeom\n";
if ($ngeom==0) {
close CON;
return;
}
print "\n";
print "For each geometry in the scan, the unconstrained degrees of freedom will be\n";
print "optimised in a set of steps. You must specify the maximum number of such\n";
print "optimisation steps allowed.\n";
print "Enter the maximum number of steps.\n";
$maxsteps=<STDIN>;
chomp $maxsteps;
open (STEP,">OPTSTEPS");
print STEP "The maximum number of optimisation steps is\n";
print STEP "$maxsteps\n";
close STEP;
print "\n";
print CON "Enter the number of bond lengths to be constrained on the path\n";
print "Enter the number of bond lengths to be constrained on the path\n";
$nbonds=<STDIN>;
chomp $nbonds;
print CON "$nbonds\n";
#if ($nbonds==0) {
#close CON;
#return;
#}
if ($nbonds >= 1) {
 print "Enter the atom numbers, initial bond length, and increment for each bond constraint\n";
 print CON "Enter the atom numbers, initial bond length, and increment for each bond constraint\n";
 for ($il=1; $il<=$nbonds; $il++) {
 $line=<STDIN>;
 chomp $line;
 print CON "$line\n";
 }
}
print CON "Enter the number of angles to be constrained\n";
print "Enter the number of angles to be constrained\n";
$nangle=<STDIN>;
chomp $nangle;
print CON "$nangle\n";
if ($nangle >= 1) {
 print "Enter the 3 atom numbers, initial angle and increment in degrees for each angle\n";
 print CON "Enter the 3 atom numbers, initial angle and increment in degrees for each angle\n";
 for ($il=1; $il<=$nangle; $il++) {
 $line=<STDIN>;
 chomp $line;
 print CON "$line\n";
 }
}
print "Enter the number of dihedral angles to be constrained\n";
print CON "Enter the number of dihedral angles to be constrained\n";
$ndiheds=<STDIN>;
chomp $ndiheds;
print CON "$ndiheds\n";
if ($ndiheds >= 1) {
 print "Enter the 4 atom numbers, initial angle and increment in degrees for each dihedral\n";
 print CON "Enter the 4 atom numbers, initial angle and increment in degrees for each dihedral\n";
 for ($il=1; $il<=$ndiheds; $il++) {
 $line=<STDIN>;
 chomp $line;
 print CON "$line\n";
 }
}
close CON;
}

sub callPrepare {
system("rm -f OUT_NONMETALCHARGES");
system "$EXEDIR/Preparegeom < molxyz > molcharge";

}

sub utilities {

    while (1) {  # always true
     &menu_init($numbered_flag, "SMFA Control ",0,
                ,"\n-submenu: Utility Programs");
     &menu_item ("Frequencies with isotopic substitution","Frequency");
     &menu_item ("Isodesmic, homodesmotic and analogous reactions","isodesmicetc");
     &menu_item ("Combining isodesmic, homodesmotic etc reactions","subtractreactions");
     &menu_item ("Electrostatic potential on the solvent-accessible-surface","solsurfpot");
     &menu_item ("Dipole Polarizability","maketotalpolar");
     &menu_item ("Dipole Hyperpolarizability","makehyperpolar_gau");
     &menu_item ("Internal Coordinates","internal");
     &menu_item ("Add H atoms","newHatoms");
     &menu_item ("Exit","Exit");

     $sel = &menu_display("",$menu_default_row,$menu_default_top,
                             $menu_default_col);

     if ($sel eq "%EMPTY%" ) { die "not enough screen lines \n";}
     if ($sel eq "Exit") {last;}

       if ($sel ne "%UP%" )
        { # call the corresponding subprogram (assuming existence)
         $ret=&$sel();
           if ($ret)
            {
            $window=&initscr();
            &menu_curses_application($window);
            } # end of: if ($ret)
        } # end of: if ($sel ne "%UP%" )
    } # while (1)
  } # end of: sub inp



 sub ltrim { my $s = shift; $s =~ s/^\s+//;       return $s };
 sub rtrim { my $s = shift; $s =~ s/\s+$//;       return $s };
 sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };


sub Frequency {
# a utility program to evaluate frequencies and spectrum
# for isotopically substituted molecules

system("clear");
system("stty sane");

open(TMP,">IN_ISOTOPES");
print "This utility evaluates the vibrational frequencies, intensities,\n";
print "and spectrum for isotopically substituted molecules\n";
print "\n";
print "The force constant matrix (contained in combinedderivs) and the\n";
print "dipole derivatives that were previously evaluated for this molecule\n";
print "in this subdirectory are used to evaluate the frequencies.\n";
print "\n";
print "You can choose to change the isotope for every atom of a chosen element,\n";
print "or you can change the isotope for particular chosen atoms (by number).\n";

print "Do you choose a particular element (Y/N)?\n";
$fans1=<STDIN>;
chomp $fans1;
$fans1=uc($fans1);
if ($fans1 eq "Y") {
print "Enter the elemental symbol and mass (in atomic mass units)\n";
$eleiso=<STDIN>;
chomp $eleiso;
print TMP "1\n";
print TMP "$eleiso\n";
}
if ($fans1 eq "N") {
my $ic=0;
my $again="Y";
print "For each atom to be substituted by an isotope, you must enter the atom number\n";
print "and the isotope mass (in atomic mass units)\n";
print "Enter the number and mass for each atom on one line, and finish with a\n";
print "RETURN (blank line):\n";
while ($again ne "N") {
$linein=<STDIN>;
if ($linein eq "\n") {
$again="N";
}else{
chomp $linein;
@two=split(/\s+/,$linein);
$natiso[$ic]=$two[0];
$masiso[$ic]=$two[1];
$ic=$ic+1;
}
}
$ic=$ic-1;
my $ig=$ic+1;
print TMP "2\n";
print TMP "$ig\n";
for (my $ie=0;$ie <=$ic;$ie++) {
print TMP "$natiso[$ie]      $masiso[$ie]\n";
}
}
my $width=5.0;
print "\n";
print "The IR spectrum will be simulated with Lorentzian peaks with a given\n";
print "full width at half maximum (FWHM). The default value of this width is 5 cm-1 \n";
print "Enter your chosen width or hit RETURN for the default:\n";
my $wid=<STDIN>;
if ($wid eq "\n") {
}else{
chomp $wid;
$width=$wid;
}
print TMP "$width\n";
close TMP;
system("$EXEDIR/Freq");
print "\n";
print "The isotope substituted frequencies are in file FREQUENCIES\n";
print "The Cartesian displacements for the normal modes are in file NORMAL_MODES\n";
print "The simulated spectrum is in file SPECTRUM\n";
system("sleep 5");
system("stty sane");
return 1
}

sub isodesmicetc {

system("stty sane");

$fract=`awk 'NR==2,NR==2 {print \$1}' IN_FRACT_LevX`;
chomp $fract;
$molname=`awk 'NR==1,NR==1 {print \$1}' xyzFILENAME`;
chomp $molname;
$natw=`awk 'NR==2,NR==2 {print \$1}' xyzFILENAME`;
chomp $natw;
$natw2=$natw+2;
system("head -$natw2 $molname > tmpfile");
system("mv tmpfile $molname");
system("cp $molname $molname.xyz");
system("babel $molname.xyz $molname.inchi");
system("rm $molname.xyz");
$molinchi=`awk 'NR==1,NR==1 {print \$1}' $molname.inchi`; 
chomp $molinchi;
system("cp seefrags seefrags.xyz");
system("babel seefrags.xyz seefrags.inchi");
system("rm seefrags.xyz");
open(TMP,"<signs.out");
$i1=1;
while (<TMP>) {
$sign[$i1]="$_";
chomp $sign[$i1];
$i1=$i1+1;
}
close TMP;
#make arrays for the frag coords
$ifrg=0;
open(SEE,"<seefrags");
while (<SEE>) {
 if($_ =~ /^\d+$/ ) {
 $ifrg=$ifrg+1;
 $start[$ifrg]=$.;
}
}
$ifrg=$ifrg+1;
$start[$ifrg]=$.+1;
close SEE;

open(TMP,"<seefrags.inchi");
$i1=1;
while (<TMP>) {
$store[$i1]="$_";
chomp $store[$i1];
$i1=$i1+1;
}
close TMP;
$nfrags=$i1-1;
for (my $i1=1;$i1 <= $nfrags;$i1++) {
$nmol[$i1]=1;
$used[$i1]=0;
}
$nmol[1]=1;
for (my $i1=1;$i1 < $nfrags;$i1++) {
 $i20=$i1+1;
 if($used[$i1] == 0) {
 for (my $i2=$i20;$i2 <= $nfrags;$i2++) {
  if($used[$i2] == 0) {
  if($store[$i1] eq $store[$i2]) {
  $nmol[$i1]=$nmol[$i1]+1;
  $used[$i2]=1;
 }
# end the used for I2
}
# end the i2 loop
}
# end used for i1
}
# end the i1 loop
}
$lastone=0;
for (my $i1=1;$i1 <= $nfrags;$i1++) {
 if($used[$i1] == 0) {
  if($sign[$i1] > 0) {$lastone=$i1};
 }
}
system("rm $molname.lhs.coords");
system("rm $molname.rhs.coords");
open(LHSI,">$molname.lhs.inchi");
open(LHSS,">$molname.lhs.coeff");
open(RHSI,">$molname.rhs.inchi");
open(RHSS,">$molname.rhs.coeff");
print LHSI "$molinchi\n";
print LHSS "1\n";
system("cp $molname $molname.lhs.coords");
open(TMP,">>OUT_SMFA");
print TMP "\n";
print TMP "                      Near isoenergetic reactions\n";
print TMP "\n";
print TMP "For the molecule contained in file $molname, an almost isoenergetic reaction\n";
print TMP "has been found, corresponding to a Level = $fract fragmentation.\n";
print TMP "\n";
print TMP "\n";
print TMP "The molecule in $molname \($molinchi\)";
print TMP "\n";
print TMP "\n";
for (my $i1=1;$i1 <= $nfrags;$i1++) {
 if (($used[$i1] == 0) && ($sign[$i1] < 0)) {
  print TMP "  +\n";
  print TMP "\n";
  print TMP " $nmol[$i1]   X   Fragment $i1    \($store[$i1]\) \n";
  print LHSI "$store[$i1]\n";
  print LHSS "$nmol[$i1]\n";
$i2=$i1+1;
$ib=$start[$i2]-1;
  system("awk ' NR=='$start[$i1]',NR == '$ib' {print  \$0} ' seefrags >> $molname.lhs.coords");
  
 }
}
  print TMP "\n";
  print TMP "  <==>\n";
  print TMP "\n";
for (my $i1=1;$i1 <= $nfrags;$i1++) {
 if (($used[$i1] == 0) && ($sign[$i1] > 0)) {
  print TMP  " $nmol[$i1]   X   Fragment $i1    \($store[$i1]\) \n";
  print RHSI "$store[$i1]\n";
  print RHSS "$nmol[$i1]\n";
$i2=$i1+1;
$ib=$start[$i2]-1;
  system("awk ' NR=='$start[$i1]',NR == '$ib' {print  \$0} ' seefrags >> $molname.rhs.coords");
 if($i1 < $lastone) {
  print TMP "  +\n";
  print TMP "\n";
 }
 }
}
close LHSI;
close LHSS;
close RHSI;
close RHSS;
system("obabel -ixyz $molname.lhs.coords -osvg -O $molname.lhs.svg");
system("obabel -ixyz $molname.rhs.coords -osvg -O $molname.rhs.svg");

print TMP "\n";
print TMP " The Cartesian coordinates for each Fragment are listed\n";
print TMP " in the file seefrags.\n";
  print TMP "\n";
print TMP " In addition, files describing the lhs and rhs of the reaction have been created\n";
print TMP " under the file names\n";
print TMP "\n";
print TMP " $molname.lhs.inichi\n";
print TMP " $molname.lhs.coeff\n";
print TMP " $molname.lhs.coords\n";
print TMP " $molname.lhs.svg\n";
print TMP " $molname.rhs.inichi\n";
print TMP " $molname.rhs.coeff\n";
print TMP " $molname.rhs.coords\n";
print TMP " $molname.rhs.svg\n";
print TMP "\n";
print TMP " Files with '.inichi' list the Inchi for each compound\n";
print TMP " Files with '.coeff' list the coefficient of the compound in the reaction formula\n";
print TMP " Files with '.coords' list the Cartesian coordinates of the compounds, visualised with VMD etc\n";
print TMP " Files with '.svg can be opended with a web browser like FireFox to give a table of\n";
print TMP " of 2D drawings of the compounds\n";
print TMP "\n";
print TMP "The $molname.lhs... and $molname.rhs...  files can be used with another utility in SMFA\n";
print TMP "to construct more 'near isoenergetic' reactions - see the User's manual\n";
close TMP;
system("clear");
print "The near isoenergetic reaction has been evaluated, and the results\n";
print "have been appended to OUT_SMFA\n";
print "\n";
print " In addition, files describing the lhs and rhs of the reaction have been created\n";
print " under the file names\n";
print "\n";
print " $molname.lhs.inichi\n";
print " $molname.lhs.coeff\n";
print " $molname.lhs.coords\n";
print " $molname.lhs.svg\n";
print " $molname.rhs.inichi\n";
print " $molname.rhs.coeff\n";
print " $molname.rhs.coords\n";
print " $molname.rhs.svg\n";
print "\n";
print " Files with '.inichi' list the Inchi for each compound\n";
print " Files with '.coeff' list the coefficient of the compound in the reaction formula\n";
print " Files with '.coords' list the Cartesian coordinates of the compounds, visualised with VMD etc\n";
print " Files with '.svg can be opended with a web browser like FireFox to give a table of\n";
print " of 2D drawings of the compounds\n";
print "\n";
print "The $molname.lhs... and $molname.rhs...  files can be used with another utility in SMFA\n";
print "to construct more 'near isoenergetic' reactions - see the User's manual\n";
system("sleep 10");
}

sub subtractreactions {
system("stty sane");

print "                       Combining near Isoenergetic reactions\n";
print "\n";
print " If you have created the files for more than one isodesmic, homodesmotic (etc) reaction,\n";
print " you can use this utility to 'subtract' one reaction from another.\n";
print "\n";
print " For example, you may have found an isodesmic, homodemotic etc reaction for each\n";
print " of two isomers. Since both reactions are near isoenergetic, a reaction that we form\n";
print " by subtracting one reaction from the other is also near isoenergetic. In this\n";
print " example, the resulting reaction could be used to estimate the energy difference\n";
print " between the two isomers. If the isomers have many similarities in structure, then\n";
print " 'subtracting' one reaction from the other will result in the cancelation of many\n";
print " components of the reactions. Thus a simpler reaction for estimating the energy\n";
print " difference betwwen the two isomers will be obtained.\n";
print "\n";
print " Many other applications for constructing near isoenergetic reactions may come to mind.\n";
print "\n";
print " The isodesmic (etc) utility created several files for the reactions with names like\n";
print " name1.lhs.inchi and name2.lhs.inchi\n";
print " where name1 and name2 were the names of the coordinate files for these molecules.\n";
print "\n";
print " All you have to do is enter these two file names, eg name1 and name2\n";
print "\n";
print " Enter the first file name\n";
$name1=<STDIN>;
chomp $name1;
$name1=~ s/^\s+|\s+$//g;
unless ( -s "$name1.lhs.coeff" ) {
print "File does not exist\n";
system("sleep 5");
return;
}
print " Enter the second file name\n";
$name2=<STDIN>;
chomp $name2;
$name2=~ s/^\s+|\s+$//g;
unless ( -s "$name2.lhs.coeff" ) {
print "File does not exist\n";
system("sleep 5");
return;
}

open(LHS,"$name1.lhs.coeff");
$i1=0;
while(<LHS>) {
$i1=$i1+1;
$lhssigns[$i1]=$_;
chomp $lhssigns[$i1];
}
close LHS;
open(LHS,"$name2.rhs.coeff");
while(<LHS>) {
$i1=$i1+1;
$lhssigns[$i1]=$_;
chomp $lhssigns[$i1];
}
close LHS;
$ntermslhs=$i1;

open(LHS,"$name1.rhs.coeff");
$i1=0;
while(<LHS>) {
$i1=$i1+1;
$rhssigns[$i1]=$_;
chomp $rhssigns[$i1];
}
close LHS;
open(LHS,"$name2.lhs.coeff");
while(<LHS>) {
$i1=$i1+1;
$rhssigns[$i1]=$_;
chomp $rhssigns[$i1];
}
close LHS;
$ntermsrhs=$i1;

open(LHS,"$name1.lhs.inchi");
$i1=0;
while(<LHS>) {
$i1=$i1+1;
$lhsinchi[$i1]=$_;
chomp $lhsinchi[$i1];
}
close LHS;
open(LHS,"$name2.rhs.inchi");
while(<LHS>) {
$i1=$i1+1;
$lhsinchi[$i1]=$_;
chomp $lhsinchi[$i1];
}
close LHS;

open(LHS,"$name1.rhs.inchi");
$i1=0;
while(<LHS>) {
$i1=$i1+1;
$rhsinchi[$i1]=$_;
chomp $rhsinchi[$i1];
}
close LHS;
open(LHS,"$name2.lhs.inchi");
while(<LHS>) {
$i1=$i1+1;
$rhsinchi[$i1]=$_;
chomp $rhsinchi[$i1];
}
close LHS;

open(LHS,"$name1.lhs.coords");
$i1=0;
while(<LHS>) {
$i1=$i1+1;
$lhscoords[$i1]=$_;
chomp $lhscoords[$i1];
}
close LHS;
open(LHS,"$name2.rhs.coords");
while(<LHS>) {
$i1=$i1+1;
$lhscoords[$i1]=$_;
chomp $lhscoords[$i1];
}
close LHS;
$ncoordslhs=$i1;

open(LHS,"$name1.rhs.coords");
$i1=0;
while(<LHS>) {
$i1=$i1+1;
$rhscoords[$i1]=$_;
chomp $rhscoords[$i1];
}
close LHS;
open(LHS,"$name2.lhs.coords");
while(<LHS>) {
$i1=$i1+1;
$rhscoords[$i1]=$_;
chomp $rhscoords[$i1];
}
close LHS;
$ncoordsrhs=$i1;

# get line numbers for frags in coords
$ic=0;
for ($i1=1;$i1 <= $ncoordslhs;$i1++) {
 if($lhscoords[$i1] =~ /^\d+$/ ) {
  $ic=$ic+1;
  $startlhs[$ic]=$i1;
 }
}
$startlhs[$ic+1]=$ncoordslhs+1;
$ic=0;
for ($i1=1;$i1 <= $ncoordsrhs;$i1++) {
 if($rhscoords[$i1] =~ /^\d+$/ ) {
  $ic=$ic+1;
  $startrhs[$ic]=$i1;
 }
}
$startrhs[$ic+1]=$ncoordsrhs+1;

# now we have both sides of the reaction, we can look for cancelations

for ($i1=1;$i1 <= $ntermslhs;$i1++) {

 for($i2=1;$i2 <= $ntermsrhs;$i2++) {

  if($lhsinchi[$i1] eq $rhsinchi[$i2]) {
   if($lhssigns[$i1] >= $rhssigns[$i2]) {
    $lhssigns[$i1]=$lhssigns[$i1]-$rhssigns[$i2];
    $rhssigns[$i2]=0;
   } else {
    $rhssigns[$i2]=$rhssigns[$i2]-$lhssigns[$i1];
    $lhssigns[$i1]=0;
   }
  }
 }
}

open(LHSS,">$name1-$name2.lhs.coeff");
open(LHSI,">$name1-$name2.lhs.inchi");
open(LHSC,">$name1-$name2.lhs.coords");
open(RHSS,">$name1-$name2.rhs.coeff");
open(RHSI,">$name1-$name2.rhs.inchi");
open(RHSC,">$name1-$name2.rhs.coords");

for ($i1=1;$i1 <= $ntermslhs;$i1++) {

if($lhssigns[$i1] > 0) {
 print LHSS "$lhssigns[$i1]\n";;
 print LHSI "$lhsinchi[$i1]\n";
 for ($j1=$startlhs[$i1];$j1 <= $startlhs[$i1+1]-1;$j1++) {
  print LHSC "$lhscoords[$j1]\n";
 }
}
}

for ($i1=1;$i1 <= $ntermsrhs;$i1++) {

if($rhssigns[$i1] > 0) {
 print RHSS "$rhssigns[$i1]\n";
 print RHSI "$rhsinchi[$i1]\n";
 for ($j1=$startrhs[$i1];$j1 <= $startrhs[$i1+1]-1;$j1++) {
  print RHSC "$rhscoords[$j1]\n";
 }
}
}

close LHSS;
close LHSI;
close LHSC;
close RHSS;
close RHSI;
close RHSC;

system("obabel -ixyz $name1-$name2.lhs.coords -osvg -O $name1-$name2.lhs.svg");
system("obabel -ixyz $name1-$name2.rhs.coords -osvg -O $name1-$name2.rhs.svg");

print "The files for the reaction formed from\n";
print "reaction $name1 minus reaction $name2\n";
print "have been written to files named $name1-$name2\n";
print "These $name1-$name2.inchi (etc) files have exactly the same structure as the original\n";
print "$name1.inchi (etc) files, and so are available for further manipulations\n";
print "\n";
print "A note of this calculations has been appended to OUT_SMFA\n";
print "The program will return to the main menu\n";
open(TMP,">>OUT_SMFA");
print TMP "The files for the reaction formed from\n";
print TMP "reaction $name1 minus reaction $name2\n";
print TMP "have been written to files named $name1-$name2\n";
print TMP "These $name1-$name2.inchi (etc) files have exactly the same structure as the original\n";
print TMP "$name1.inchi (etc) files, and so are available for further manipulations\n";
print TMP "\n";
close TMP;

system("sleep 10");

}





sub maketotalpolar {
system("clear");
system("stty sane");

$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
$n++;
}
open(TMP,"<IN_JOBTYPE");
while (<TMP>) {
if ($. == 2) {$deriv="$_"};
}
chomp $deriv;

$qmprog=$store[1];

# first use the group polarizabilities

if (-s "Polarisabilities.out") {
for (my $ic=1;$ic <= 7;$ic++) {
$totpolar[$ic]=0;
}

open(TMP,"<Polarisabilities.out");
while (<TMP>) {
$line1="$_";
@line=split(/\s+/,$line1);
for (my $ic=1;$ic <= 7;$ic++) {
 $totpolar[$ic]=$totpolar[$ic]+$line[$ic];
}
}
close TMP;
open(TMP,">>OUT_SMFA");
print TMP "\n";
print TMP "The total molecular dipole polarizability is approximated\n";
print TMP "by a sum over group polarizabilities as:\n";
print TMP "\n";
print TMP "     xx         yx         yy         zx         zy         zz\n";
print TMP "\n";
printf TMP "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",$totpolar[1],$totpolar[2],$totpolar[3],$totpolar[4],$totpolar[5
],$totpolar[6];
print TMP "\n";
printf TMP "%s   %10.3f\n",'The isotropic polarizability = ',$totpolar[7];
print TMP "\n";
}


if (($qmprog == 1) || ($qmprog == 3)) {
print "SMFA can calculate the dipole polarizability as a sum over the fragment polarizabilities,\n";
print "only for GAUSSIAN or QCHem.\n";
print TMP "SMFA can calculate the dipole polarizability as a sum over the fragment polarizabilities,\n";
print TMP "only for GAUSSIAN or QCHem.\n";
close TMP;
system("sleep 5");
return;
}

my $nfrgs = `ls COORD* | wc -l`;
chomp $nfrgs;
open(SG,"<signs.out");
while(<SG>) {
$sign[$.]=$_;
chomp $sign[$.];
$sign[$.]=ltrim($sign[$.]);
}
close SG;
for (my $ic=1;$ic <= 10;$ic++) {
$totpolar[$ic]=0;
}

if ( $qmprog == 2 ) {

$pol=`grep 'POLAR' FRAG1.com | wc -l`;
chomp($pol);
$pol=~ s/^\s+|\s+$//g;if ( $pol == 0 && $deriv != 2 ) {
print "There was no POLAR or FREQ calculation performed. Returning to the menu.\n";
system("sleep 5");
return;
}
print TMP "For GAUSSIAN, SMFA can calculate the dipole polarizability as a sum over the fragment polarizabilities.\n";
print TMP "This is a more accurate approximation.\n";

for (my $ic=1;$ic <= $nfrgs;$ic++) {
&getpolar_gau("FRAG$ic.log");
for (my $ie=1;$ie <= 6;$ie++) {
 $totpolar[$ie]=$totpolar[$ie]+$sign[$ic]*$Polar[$ie];
}
$totpolar[7]=$totpolar[7]+$sign[$ic]*($Polar[1]+$Polar[3]+$Polar[6]);
}
$totpolar[7]=$totpolar[7]/3;
print TMP "\n";
print TMP "     xx         yx         yy         zx         zy         zz\n";
printf TMP "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",$totpolar[1],$totpolar[2],$totpolar[3],$totpolar[4],$totpolar[5],$totpolar[6];
print TMP "\n";
printf TMP "%s   %10.3f\n",'The isotropic polarizability = ',$totpolar[7];
}

if ( $qmprog == 4) {

$pol=`grep 'MOPROP' IN_QCH | wc -l`;
chomp($pol);
$pol=~ s/^\s+|\s+$//g;
if ( $pol == 0 ) {
print "The \$rem entry 'MOPROP' was not included. Returning to the menu.\n";
system("sleep 5");
return;
}
print TMP "For QChem, SMFA can calculate the dipole polarizability as a sum over the fragment polarizabilities.\n";
print TMP "This is a more accurate approximation.\n";
for (my $ic=1;$ic <= $nfrgs;$ic++) {
&getpolar_qch("FRAG$ic.log");
for (my $ie=1;$ie <= 9;$ie++) {
 $totpolar[$ie]=$totpolar[$ie]+$sign[$ic]*$Polar[$ie];
}
$totpolar[10]=$totpolar[10]+$sign[$ic]*($Polar[1]+$Polar[5]+$Polar[9]);
}
$totpolar[10]=$totpolar[10]/3;
print TMP "\n";
print TMP "     xx         yx         yy         zx         zy         zz\n";
# reverse an apparent sign error in QChem
#for (my $ic=1;$ic <= 10;$ic++) {
#$totpolar[$ic]=-$totpolar[$ic];
#}
printf TMP "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",$totpolar[1],$totpolar[4],$totpolar[5],$totpolar[7],$totpolar[8]
,$totpolar[9];
print TMP "\n";
printf TMP "%s   %10.3f\n",'The isotropic polarizability = ',$totpolar[10];
print TMP "\n";
}



close TMP;

print "The polarizability has been evaluated and the result has been appended to OUT_SMFA\n";
system("sleep 4");

}

sub solsurfpot {
system("clear");
system("stty sane");

print "             Electrostatic potential on the solvent-accessible-surface\n";
print "\n";
print "\n";
print " SMFA will generate the data to allow you to visualise the electrostatic potential\n";
print " on the solvent accessible surface surrounding the molecule.\n";
print " The electrostatic potential is generated from a charge and dipole, quadrupole,\n";
print " octapole and hexadecapole moments distributed onto every atom in the molecule.\n";
print " These charges and electrostatic moments are generated from the electron densities that\n";
print " were generated for the Level = 1 fragments, using Stone's method applied to the output\n";
print " from GAUSSIAN or QChem. Similarly, Stone's method is integrated into GAMESS, and the\n";
print " corresponding moments are used. For NWChem, only distributed charges on the atoms are\n";
print " available for use.\n";
print "\n";
print " The solvent accessible surface is defined as the surface that is separated from the\n";
print " nearest atom in the molecule by the Van der Waals radius of that atom plus a radius\n";
print " associated with the solvent molecule.\n";
print "\n";
print " This utility program provides two output files for graphics:\n";
print " The first produces a file called SOLACCSURFACE. This is an 'xyz' format file\n";
print " which can be loaded by programs like VMD or MacMolPlt to show you the\n";
print " solvent-accessible-surface surrounding the molecule. This graphic is useful in\n";
print " conjunction with the second file, called ESP.cube. This file has the format of a\n";
print " GAUSSIAN 'cubegen' file, and can be read by VMD (and other programs) to show\n";
print " the electrostatic potential on the solvent-accessible-surface. See the User's manual\n";
print " for a guide to using VMD with this file.\n";
print "\n";
print "You must enter a value for this solvent radius (the default is 1.4 Angstrom)\n";
$vdwrad=<STDIN>;
if($vdwrad eq "\n") {
$vdwrad=1.4;
}else{
chomp $vdwrad;
}
print "You must enter the desired density of points on the surface (per square Angstrom)\n";
print " The default value is 1.0\n";
$rho=<STDIN>;
if($rho eq "\n") {
$rho=1.0;
}else{
chomp $rho;
}
print " You must enter the number of grid points for the 'cube' of ESP data,\n";
print " The recommended value is 100 (per axis direction)\n";
$npoints=<STDIN>;
chomp $npoints;
print " A note on these calculations will be appended to OUT_SMFA\n";
open(TMP,">IN_SOLSURF");
print TMP "$vdwrad\n";
print TMP "$rho\n";
print TMP "$npoints\n";
close TMP;
open(TMP,">>OUT_SMFA");
print TMP "\n";
print TMP "\n";
print TMP "\n";
print TMP " SMFA will generate the data to allow you to visualise the electrostatic potential\n";
print TMP " on the solvent accessible surface surrounding the molecule.\n";
print TMP " The electrostatic potential is generated from a charge and dipole, quadrupole,\n";
print TMP " octapole and hexadecapole moments distributed onto every atom in the molecule.\n";
print TMP " These charges and electrostatic moments are generated from the electron densities that\n";
print TMP " were generated for the Level = 1 fragments, using Stone's method applied to the output\n";
print TMP " from GAUSSIAN or QChem. Similarly, Stone's method is integrated into GAMESS, and the\n";
print TMP " corresponding moments are used. For NWChem, only distributed charges on the atoms are\n";
print TMP " available for use.\n";
print TMP "\n";
print TMP " The solvent accessible surface is defined as the surface that is separated from the\n";
print TMP " nearest atom in the molecule by the Van der Waals radius of that atom plus a radius\n";
print TMP " associated with the solvent molecule.\n";
print TMP "\n";
print TMP " This utility program provides two output files for graphics:\n";
print TMP " The first produces a file called SOLACCSURFACE. This is an xyz format file\n";
print TMP " which can be loaded by programs like VMD or MacMolPlt to show you the\n";
print TMP " solvent-accessible-surface surrounding the molecule. This graphic is useful in\n";
print TMP " conjunction with the second file, called ESP.cube. This file has the format of a\n";
print TMP " GAUSSIAN cubegen file, and can be read by VMD (and other programs) to show\n";
print TMP " the electrostatic potential on the solvent-accessible-surface. See the User's manual\n";
print TMP " for a guide to using VMD with this file.\n";
print TMP "\n";

print TMP " The solvent radius entered is $vdwrad\n";
print TMP "\n";
print TMP " The requested density (per square Angstrom) of points on the surface is $rho\n";
print TMP "\n";
print TMP "The number of grid points for the 'cube' of ESP data is $npoints per axis\n";
print TMP "\n";
close TMP;
system("$EXEDIR/AllMoments");
system("$EXEDIR/SolAccSurface");
system("$EXEDIR/makecube");

}

sub viewfrags {

open(SEE,">seefrags");

my $nfrgs = `ls COORD* | wc -l`;
for (my $ic=1; $ic<=$nfrgs; $ic++) {
my $mf=`wc -l COORD$ic`;
my @fields = split / /, $mf;
$natomf=$fields[0];
$natomf = $natomf - 1;
print SEE "$natomf\n";
print SEE "Fragment $ic\n";
close SEE;
system("tail -$natomf COORD$ic >> seefrags");
open(SEE,">>seefrags");
}
close SEE;
}


sub optimise {

# tidy up
system("rm Newco*");
system("rm hessian*");
system("rm gradient*");
system("rm energy*");
system("rm combine*");
system("rm -f CONVERGEDCOORDS");
system("rm -f UNCONVERGEDCOORDS");

# get some parameters
$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while (<$fh>){
if(/The job type is/) {$nsave=$n+1};
$row="$_";
chomp $row;
$store[$n]=$row;
$n++;
}
$n=$n-1;
close $fh;
$nproc0=$store[$n];
for (my $ic=0; $ic<=$n; $ic++) {
 if ("$store[$ic]" eq "Is long range dispersion accounted for?")  {
  $ie=$ic+1;
  $disp=$store[$ie];
 }
}
$deriv_orig=$store[$nsave];

$nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;

$maxsteps=`awk 'NR==2,NR==2 {print \$1}' OPTSTEPS`;

#here because $deriv_orig=3 (opt)

$package =`awk 'NR==1,NR==1 {print \$0}' IN_PACKAGE`;
chomp $package;
$qmprog=$package;;

# first the hessian

system("cp IN_JOBTYPE IN_JOBTYPE_original");

open(TMP,">IN_JOBTYPE");
print TMP "jobtype file for opt\n";
print TMP "2\n";
close TMP;
$deriv=2;
open(TMP,">STEPNUMBER");
print TMP "0\n";
close TMP;


if($package ==1) {
# old bipassed code
#&SMFAgaminputs_MAC;
#&SMFArungam;
#&dodaltonpolar;
#&dodaltondisp;
#&extract_gam;
#&makenbgdma_gam;
#&ReadPolar_dal;
}

if($package ==2) {
system("cp IN_G09 IN_G09_original");
system("cp IN_G09_step0 IN_G09");
&Lev0_chg_MAC_iter_gau;
&SMFAgauinputs_MAC;
&dodaltondisp;
&SMFArunallgau;
&anal_MAC;
}

if($package ==3) {
# old bipassed code
#system("cp NWCcommands NWCcommands_original");
#system("cp NWCcommands tmpfile");
#system("sed 's/OPTIMIZE/GRADIENT/' tmpfile > NWCcommands");
#system("echo 'TASK SCF HESSIAN' >> NWCcommands");
#&Lev0_chg_MAC_iter_nwc;
#&SMFAnwcinputs_MAC;
#&SMFArunnwc;
#&dodaltonpolar;
#&dodaltondisp;
#&extract_nwc;
#&makenbgdma_nwc;
#if ("$disp" eq "Y") {
# &ReadPolar_dal;
# } else {
#   if ($nchgs eq 0) {
#    &ReadPolar_dal;
#   }
#}
}


if($package ==4) {
system("cp IN_QCH IN_QCH_original");
system("cp IN_QCH_step0 IN_QCH");
&Lev0_chg_MAC_iter_qch;
&SMFAqchinputs_MAC;
&dodaltondisp;
&SMFArunallqch;
&anal_MAC;
}

system("cp combinedderivs combinedderivs.0");

# now have the gradient and hessian for the first step
# store the energy,gradient and hessian
&getcombdata("0");

$natom=`awk 'NR==1,NR==1 {print \$1}' name.xyz`;
chomp $natom;
$n3= 3 * $natom;
$rad= 0.005 * $n3;
open(TMP,">TRUSTRADIUS");
print TMP "$rad\n";
close TMP;


system("rm optout");
system("$EXEDIR/ROGEROPT >> optout");
system("mv deltax deltax.0");
system("cp fort.67 fort.67_0");


system("cp xyzFILENAME xyzFILENAME_original");
$numberofatoms=`awk 'NR==2,NR==2 {print \$1}' xyzFILENAME`;
chomp $numberofatoms;
open(TMP,">xyzFILENAME");
print TMP "Newcoords.xyz\n";
print TMP "$numberofatoms\n";
close TMP;

# return deriv to 1 (gradient only)
open(TMP,">IN_JOBTYPE");
print TMP "jobtype file for opt\n";
print TMP "1\n";
close TMP;
$deriv=1;

# update the fragments
system("rm -f WARNINGS");
system("$EXEDIR/Preparegeom < Newcoords.xyz > rub");
system("cat -s WARNINGS >> optout");
system("$EXEDIR/SOLCH");
open(TMP2,">READY2GO");
print TMP2 "1\n";
close(TMP2);

&fragonly;


if($package == 1) {

}

if($package == 2) {
system("cp IN_G09_original tmpfile");
system("sed 's/OPT/FORCE/' tmpfile > IN_G09");
}
if($package == 3) {
system("cp NWCcommands_original NWCcommands");
system("cp NWCcommands tmpfile");
system("sed 's/OPTIMIZE/GRADIENT/' tmpfile > NWCcommands");
}
 
if($package == 4) {
system("cp IN_QCH_original tmpfile");
system("sed 's/JOBTYPE OPT/JOBTYPE FORCE/' tmpfile > IN_QCH");
}



# now take the steps

for ($istep=1; $istep < $maxsteps; $istep++) {

open(TMP,">STEPNUMBER");
print TMP "$istep\n";
close TMP;

if($package ==1) {
#&SMFAgaminputs_MAC;
#&SMFArungam;
#&dodaltonpolar;
#&dodaltondisp;
#&extract_gam;
#&makenbgdma_gam;
#&ReadPolar_dal;
}
if($package ==2) {
&Lev0_chg_MAC_iter_gau;
&SMFAgauinputs_MAC;
&dodaltondisp;
&SMFArunallgau;
&anal_MAC;
}
if($package ==3) {
#&Lev0_chg_MAC_iter_nwc;
#&SMFAnwcinputs_MAC;
#&SMFArunnwc;
#&dodaltonpolar;
#&dodaltondisp;
#&extract_nwc;
#&makenbgdma_nwc;
#if ("$disp" eq "Y") {
# &ReadPolar_dal;
# } else {
#   if ($nchgs eq 0) {
#    &ReadPolar_dal;
#   }
#}
}
if($package ==4) {
&Lev0_chg_MAC_iter_qch;
&SMFAqchinputs_MAC;
&dodaltondisp;
&SMFArunallqch;
&anal_MAC;
}

# store the energy, gradient and hessian

# update the hessian
system("$EXEDIR/updatehessian");
&getcombdata("$istep");
system("cp combinedderivs combinedderivs.$istep");

system("$EXEDIR/ROGEROPT >> optout");
system("mv deltax deltax.$istep");
system("cp Newcoords.xyz Newcoords.xyz.$istep");
system("cp fort.67 fort.67_$istep");

$converged=`grep converged optout | wc -l`;
chomp $converged;
if($converged ne "0") {
system("cp Newcoords.xyz CONVERGEDCOORDS");
open(TMP,">>OUT_SMFA");
print TMP "\n";
print TMP "The geometry optimsation has converged in $istep steps.\n";
print TMP "The converged geometry is contained in the file CONVERGEDCOORDS.\n";
print TMP "Intermediate output from the optimisation process can be viewed\n";
print TMP "in the file optout";
print TMP "\n";
print TMP "The initial geometry had a total energy of\n";
$energy0=`awk 'NR==1,NR==1 {print \$1}' energy.0`;
chomp $energy0;
print TMP "$energy0\n";
print TMP "\n";
print TMP "The final geometry has a total energy of\n";
$istep1=$istep-1;
$energylast=`awk 'NR==1,NR==1 {print \$1}' energy.$istep1`;
chomp $energylast;
print TMP "$energylast\n";
close TMP;
system("cp IN_JOBTYPE_original IN_JOBTYPE");
system("cp xyzFILENAME_original xyzFILENAME");
if($package == 2) {
system("cp IN_G09_original IN_G09");
}
if($package == 3) {
system("cp NWCcommands_original NWCcommands");
}
if($package == 4) {
system("cp IN_QCH_original IN_QCH");
}
return;
}


# update the fragments
system("$EXEDIR/Preparegeom < Newcoords.xyz > rub");
system("cat -s WARNINGS >> optout");
system("$EXEDIR/SOLCH");
open(TMP2,">READY2GO");
print TMP2 "1\n";
close(TMP2);

&fragonly;

open(TMP,">>OUT_SMFA");
print TMP "\n";
print TMP "Intermediate geometry $istep during the optimisation\n";
print TMP "\n";
close TMP;

$maxsteps1=$maxsteps-1;
if ($istep==$maxsteps1) {
system("cp Newcoords.xyz UNCONVERGEDCOORDS");
open(TMP,">>OUT_SMFA");
print TMP "The optimisation did not converge\n";
print TMP "The final geometry is contained in the file UNCONVERGEDCOORDS\n";
print TMP "The final geometry has a total energy of - unconverged\n";
$energylast=`awk 'NR==1,NR==1 {print \$1}' energy.$istep`;
chomp $energylast;
print TMP "$energylast\n";
close TMP;
if($package == 2) {
system("cp IN_G09_original IN_G09");
}
if($package == 3) {
system("cp NWCcommands_original NWCcommands");
}
if($package == 4) {
system("cp IN_QCH_original IN_QCH");
}
return;
}

# end the loop over steps
}

# end optimise
}


sub getcombdata {
$flag="$_[0]";
open(TMP,"<combinedderivs");
$skip0=0;
$writit0=0;
$skip1=0;
$writit1=0;
$skip2=0;
$writit2=0;
open(ENER,">energy.$flag");
open(GRAD,">gradient.$flag");
open(HESS,">hessian.$flag");
while(<TMP>) {
if(/The energy is/) {$skip0 = $. + 1};
if($skip0 == $.) {$writit0 = 1};
if(/The coordinates are/) {$writit0=0};
if($writit0 == 1) {
$energy="$_";
chomp $energy;
print ENER "$energy\n";
}
if(/First derivatives/) {
$skip1 = $. + 1;
}
if($skip1 == $.) {$writit1 = 1};
if(/Upper triangle of the second derivatives/) {
$skip2 = $. + 1;
$writit1=0;
}
if($writit1 == 1) {
$grad="$_";
chomp $grad;
print GRAD "$grad\n";
}
if($skip2 == $.) {$writit2 = 1};
if($writit2 == 1) {
$hess="$_";
chomp $hess;
print HESS "$hess\n";
}
}
close ENER;
close GRAD;
close HESS;
}

sub scan {
# the subroutine for a scan
system("stty sane");

system("cp xyzFILENAME xyzFILENAME_scanstart");
system("cp IN_JOBTYPE IN_JOBTYPE_scanstart");

$nscan=`awk 'NR==2,NR==2 {print \$1}' IN_OPTSCAN`;
chomp $nscan;
if ($nscan==0) {
return;
}

system("rm -f IN_CONSTRAINTS*");
system("$EXEDIR/makeIN_CONSTRAINTS");

open(TMP,">>OUT_SMFA");
print TMP "\n";
print TMP " Performing an optimised scan\n";
print TMP "\n";
close TMP;

system("rm -f optout*");
system("rm -f SCANcon*");
system("rm -f SCANuncon*");

for ($isc=1; $isc <= $nscan; $isc++) {

open(TMP1,">IN_JOBTYPE");
print TMP1 "jobtype file for opt\n";
print TMP1 "2\n";
close TMP1;
$deriv=2;
system("cp IN_CONSTRAINTS_$isc IN_CONSTRAINTS");
open(TMP,">>OUT_SMFA");
print TMP "Scan geometry $isc with the following constraints\n";
close TMP;
system("cat IN_CONSTRAINTS >> OUT_SMFA");

&optimise;

$converged=`grep converged optout | wc -l`;
chomp $converged;
if($converged ne "0") {
system("cp Newcoords.xyz SCANcon$isc");
$numberofatoms=`awk 'NR==2,NR==2 {print \$1}' xyzFILENAME`;
chomp $numberofatoms;
open(TMP2,">xyzFILENAME");
print TMP2 "SCANcon$isc\n";
print TMP2 "$numberofatoms\n";
close TMP2;
open(TMP,">>OUT_SMFA");
print TMP "\n";
print TMP "Scan geometry $isc converged. See SCANcon$isc and optout files";
print TMP "\n";
close TMP;
}else{
system("cp Newcoords.xyz SCANuncon$isc");
$numberofatoms=`awk 'NR==2,NR==2 {print \$1}' xyzFILENAME`;
chomp $numberofatoms;
open(TMP2,">xyzFILENAME");
print TMP2 "SCANuncon$isc\n";
print TMP2 "$numberofatoms\n";
close TMP2;
open(TMP,">>OUT_SMFA");
print TMP "\n";
print TMP "Scan geometry $isc did not converge. See SCANuncon$isc and optout files";
print TMP "\n";
close TMP;
}

system("$EXEDIR/Preparegeom < Newcoords.xyz > rub");
system("cat -s WARNINGS >> optout");
system("cp optout optout.$isc");

open(TMP2,">READY2GO");
print TMP2 "1\n";
close(TMP2);

&fragonly;

# end isc loop
}

system("cp xyzFILENAME_scanstart xyzFILENAME");
system("cp IN_JOBTYPE_scanstart IN_JOBTYPE");

# end scan
}

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
system("cat -s junk >> OUT_SMFA");
}

sub makepbsfile {
system("stty sane");
system("clear");
system("stty echo");

print " A file called pbsfile must contain the instructions\n";
print " for the PBS that controls job queues, memory requirements,\n";
print " time limits, the number of cpus requested, etc.\n";
print " This file also loads the quantum chemistry program packages.\n";
print "\n";
$localfile="pbsfile";
if ( -s $localfile ) {
 open(PBS,'<:encoding(UTF-8)',$localfile);
 $n=0;
 while (<PBS>) {
  $n=$n+1;
  $row="$_";
  chomp $row;
  $store[$n]=rtrim($row);
 }
 print " The current file of PBS instructions contains the following:\n";
 print "\n";
 for ($i=1;$i <= $n;$i++) {
 print "$store[$i]\n";
 }
 close PBS;
$tmpfile="tmppbs";
}else{
$filename="$CODEDIR/pbsfile_standard";
$tmpfile="tmppbs";
system("cp -f $filename $tmpfile");
 open(PBS,'<:encoding(UTF-8)',$tmpfile);

 $n=0;
 while (<PBS>) {
  $n=$n+1;
  $row="$_";
  chomp $row;
  $store[$n]=rtrim($row);
 }
 print " The standard file of PBS instructions contains the following:\n";
 print " NOTE that the standard file refers to a fictitious project xxx\n";
 print "\n";
 for ($i=1;$i <= $n;$i++) {
 print "$store[$i]\n";
 }
 close PBS;
}

for ($i=1;$i < 200; $i++) {
 print "\n";
 print "To add a new line, Enter A\n";
 print "To change any line in the PBS instructions, Enter C\n";
 print "To remove any line in the PBS instructions, Enter R\n";
 print "To terminate editing this file, simply hit the RETURN key\n";
 $dem=<STDIN>;
 if ($dem eq "\n") {
  last;
 }
 chomp $dem;
 $dem=uc($dem);
 if ($dem eq "A") {
  print "Enter the new line\n";
  $line=<STDIN>;
  chomp $line;
  $n=$n+1;
  $store[$n]=rtrim($line);
 }
 if ($dem eq "R") {
  print " To remove a current entry, simply type in that line (exactly)\n";
  $line=<STDIN>;
  chomp $line;
  $line=rtrim($line);
  for ($j=1;$j <= $n; $j++) {
   if ("$line" eq "$store[$j]") {
    $store[$j]="";
   }
  }
 }
 if ($dem eq "C") {
  print " First enter the line (text) you want to change\n";
  $oldline=<STDIN>;
  chomp $oldline;
  $oldline=rtrim($oldline);
  print " Then enter the corrected line (text)\n";
  $row=<STDIN>;
  chomp $row;
  for ($j=1;$j <= $n; $j++) {
   if ("$oldline" eq "$store[$j]") {
    $store[$j]=$row;
   }
  }
 }
 redo
}
 print "\n";
 print "The pbsfile now contains:\n";
 print "\n";
open (PBS,">$tmpfile");
for ($i=1;$i <= $n;$i++) {
 print PBS "$store[$i]\n";
 print "$store[$i]\n";
}
close PBS;
system("cp $tmpfile $localfile");
system("rm $tmpfile");
system("sleep 5");

#system("clear");
#$localfile="qsubfile";
#if ( -s $localfile ) {
#print "The current command for submitting SMFA for calculations is\n";
#print "\n";
#system("cat qsubfile");
#print "\n";
#print "If this is OK, simply hit RETURN, else enter a new command\n";
#  $line=<STDIN>;
#  chomp $line;
# if($line eq "") {
#  return;
# }
#system("echo $line > qsubfile");
#return;
#}

#print "You must enter the command that submits SMFA calculations to a queue,\n";
#print "for example,\n";
#print "\n";
#print "qsub -q normal\n";
#print "\n";
#print "Enter the appropriate command\n";
#  $line=<STDIN>;
#  chomp $line;
#system("echo $line > qsubfile");

}

sub internal {
system("stty sane");
system("clear");
system("stty echo");
print " To find a bond length, enter the 2 atom numbers\n";
print " To find a bond angle, enter the 3 atom numbers\n";
print " To find a dihedral angle, enter the 4 atom numbers\n";
print " To exit, hit RETURN\n";
for ($i=1;$i < 2000; $i++) {
 $line=<STDIN>;
 if($line eq "\n") {
  last;
 }
 chomp $line;
 $line=trim($line);
 open(RUB,">rub");
 print RUB "$line/n";
 close RUB;
 @num=split(/\s+/,$line);
 $length= scalar(@num);
if ($length == 2) {
system("$EXEDIR/bond < rub");
}
if ($length == 3) {
system("$EXEDIR/angle < rub");
}
if ($length == 4) {
system("$EXEDIR/dihedral < rub");
}

 redo
}
system("rm rub");
}

sub newHatoms {
system("stty sane");
system("clear");
system("stty echo");
print " To get coordinates for a H atom to replace one that is missing,\n";
print " you will need to enter the atom numbers for three atoms in the structure:\n";
print " The H atom will be bonded by 1 Angstrom to the second atom;\n";
print " the first atom is taken to lie on an imaginary +x axis from the second atom;\n";
print " the third atom is taken to lie in an imaginary xy plane,\n";
print " in the +y direction from the second atom.\n";
print " You will then enter two spherical polar angles, theta and psi\n";
print " (your best guess), and the program will output the coordinates,\n";
print " which you can paste into a modified coordinate file.\n";
for ($i=1;$i < 2000; $i++) {
print "\n";
print " Enter the 3 atom numbers, or enter RETURN to finish\n";
$line=<STDIN>;
 if($line eq "\n") {
  return;
 }
 chomp $line;
# $line=trim($line);
 open(RUB,">rub");
 print RUB "$line\n";
print " Enter the two angles in degrees\n";
$line=<STDIN>;
 if($line eq "\n") {
  return;
 }
 chomp $line;
# $line=trim($line);
 print RUB "$line\n";
 close RUB;
system("$EXEDIR/newHatom < rub");

redo
}
}


sub catELECTS {

$already=`grep 'Total numbers of atoms and electrons' OUT_SMFA  | wc -l`;
chomp $already;
if ( $already == 0 ) {
system("cat OUT_ELECTRONS_SUMMARY >> OUT_SMFA");
&largest;
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
 $Ofrags=trim($line);
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

sub restart {
system("clear");
system("stty sane");

print "You can restart the SEQUENTIAL or PARALLEL quantum chemistry calculations in SMFA\n";
print "if you have suffered a hardware failure, or similar catastrophe.\n";
print "\n";
print "Before you restart, you can change the time, memory or disc space requests.\n";
print "\n";
print "For SEQUENTIAL calculations, you can change the time, memory or disc space requests,\n";
print "by editing the pbsfile directly, or via main menu item 4, and then come back here\n";
print "to perform the restart.\n";
print "\n";
print "For PARALLEL calculations, you can change the time, memory or disc space requests,\n";
print "BUT NOT THE NUMBER OF PROCESSORS, by editing the pbsfile directly, or via main menu\n";
print "item 4, and then come back here to perform the restart.\n";
print "\n";
print "You CANNOT restart the SEQUENTIAL or PARALLEL quantum chemistry calculations in SMFA\n";
print "if you have changed anything related to the molecular structure or bonding,\n";
print "or the fragmentation parameters, or the electronic structure calculation.\n";
print "\n";
print "You CANNOT restart a GEOMETRY OPTIMISATION or SCAN using this facility\n";
print "You should store any optimised energies and geometries, and re-run the\n";
print "optimisation or scan using the last 'Newcoords.xyz' geometry as the\n";
print "starting geometry.\n";
print "\n";
print "To exit this menu item, hit RETURN\n";
print "\n";
print "To continue, enter RESTART\n";
$cont=<STDIN>;
if ( $cont eq "\n" ) {
print "RESTART ABORTED\n";
system("sleep 5");
return;
}
chomp $cont;
$cont=uc($cont);
if ( $cont ne "RESTART" ) {
print "RESTART ABORTED\n";
system("sleep 5");
return;
}
print "\n";
print "If the calculation was running SEQUENTIALLY, you must continue that way.\n";
print "If the calculation was running in PARALLEL, you must continue that way.\n";
print "Enter S or P to continue as SEQUENTIAL or PARALLEL\n";
$way=<STDIN>;
chomp($way);
$way=~ s/^\s+|\s+$//g;
$way=uc($way);
if ( $way ne "S" && $way ne "P" ) {
print "RESTART ABORTED\n";
system("sleep 5");
return;
}
if ( $way eq "P" ) {
&restartpar;
}


$out="OUTLIST";
if ( -s $out ) {
$nf=`wc -l $out`;
($nfin)=split(/\s+/,$nf);

$storelogs="storelogs";
system("mkdir $storelogs");
open(OUT,"<OUTLIST");
while (<OUT>) {
$file="$_";
chomp($file);
system("mv $file $storelogs");
}
close OUT;
}else{
print "\n";
print "There were no completed electronic structure calculations which could be saved.\n";
print "So, you should simply re-run the standard SMFA process.\n";
$cont="dead";
system("sleep 10");
return;
}

$out="FCHKLIST";
if( -s $out) {
$storefchks="storefchks";
system("mkdir $storefchks");
open(OUT,"<FCHKLIST");
while (<OUT>) {
$file="$_";
chomp($file);
system("mv $file $storefchks");
}
close OUT;
}

# now dalton files
 $out="OUTLISTDAL";
 if ( -s $out ) {
 $nf=`wc -l $out`;
 ($nfin)=split(/\s+/,$nf);

 $storelogs="storelogsdal";
 system("mkdir $storelogs");
 open(OUT,"<OUTLISTDAL");
 while (<OUT>) {
 $file="$_";
 chomp($file);
 system("mv $file $storelogs");
 }
 close OUT;
}

&viewinp;
&fragonly;
#system("rm -f INLIST");
#system("rm -f INLISTDAL");
&runall;


# end sub
}

sub createrunscripts {


#$line=`awk /ncpus=/ pbsfile`;
&getncpus;
# getncpus creates $line;

@bline=split("=",$line);
$ncpus=$bline[1];
chomp($ncpus);
$ncpus=~ s/^\s+|\s+$//g;


$manynodes=`awk 'NR==1,NR==1 {print \$1}' MULTINODE`;
chomp($manynodes);
if($manynodes eq "Y") {
 $multinode=`awk 'NR==1,NR==1 {print \$0}' MULTINODEDATA`;
 chomp($multinode);
 $multinode=~ s/^\s+|\s+$//g;
}else{
 $multinode="";
}

system("echo $ncpus >> JOBSDATA");
system("$EXEDIR/allocinput < JOBSDATA");

if($manynodes eq "Y") {
$package =`awk 'NR==1,NR==1 {print \$0}' IN_PACKAGE`;
chomp $package;
if ($package == 1) {
&multnodeinput("gaminput");
}
if ($package == 2) {
&multnodeinput("gauinput");
if ( -s "INLISTCHG" ) {
&multnodeinput("gauchinput");
}
}
if ($package == 3) {
&multnodeinput("nwcinput");
if ( -s "INLISTCHG" ) {
&multnodeinput("nwcchinput");
}
}
if ($package == 4) {
&multnodeinput("qchinput");
if ( -s "INLISTCHG" ) {
&multnodeinput("qchchinput");
}
}
}

$ncpusch=`awk 'NR==2,NR==2 {print \$1}' NPROCS`;
chomp($ncpusch);
$ncpusch=~ s/^\s+|\s+$//g;
$ncpusdal=`awk 'NR==3,NR==3 {print \$1}' NPROCS`;
chomp($ncpusdal);
$ncpusdal=~ s/^\s+|\s+$//g;

if ($ncpusch > $ncpus) {$ncpusch=$ncpus};
if ($ncpusdal > $ncpus) {$ncpusdal=$ncpus};
system("echo $ncpusch > NCPUSCH");

system("rm -f runchg");
open(RUN,">>runchg");

if ( -s "INLISTCHG" ) {
open(LIST,"<INLISTCHG");
print RUN "#!/bin/bash\n";
while (<LIST>) {
if ($. <= $ncpusch) {
$job="$_";
chomp($job);
$job=~ s/^\s+|\s+$//g;
print RUN "$CODEDIR/runlist.pl $job INLISTCHG $. &\n";
}
}
print RUN "wait\n";
close LIST;
close RUN;
$del=join("","1,",$ncpusch,"d");
system("sed -i '$del' INLISTCHG");
system("chmod +x runchg");

if($manynodes eq "Y") {
&multnodeinput("runchg");
system("chmod +x runchg");
}
}

system("rm -f runpar");
open(RUN,">>runpar");
print RUN "#!/bin/bash\n";

open(LIST,"<INLIST");
while (<LIST>) {
if ($. <= $ncpus) {
$job="$_";
chomp($job);
$job=~ s/^\s+|\s+$//g;

print RUN "$CODEDIR/runlist.pl $job INLIST $. &\n";
}
}
print RUN "wait\n";
close LIST;
close RUN;
$del=join("","1,",$ncpus,"d");
system("sed -i '$del' INLIST");
system("chmod +x runpar");
if($manynodes eq "Y") {
&multnodeinput("runpar");
system("chmod +x runpar");
}

if ( -s "INLISTDAL" ) {
system("rm -f rundal");
open(RUN,">>rundal");
print RUN "#!/bin/bash\n";
open(LIST,"<INLISTDAL");
while (<LIST>) {
if ($. <= $ncpusdal) {
$job="$_";
chomp($job);
$job=~ s/^\s+|\s+$//g;
print RUN "$CODEDIR/runlistdal.pl $job INLISTDAL $. &\n";
}
}
print RUN "wait\n";
close LIST;
close RUN;
$del=join("","1,",$ncpusdal,"d");
system("sed -i '$del' INLISTDAL");
system("chmod +x rundal");
if($manynodes eq "Y") {
&multnodeinput("rundal");
system("chmod +x rundal");
}

}

}

sub multnodeinput {
$fin="@_";

$manynodes=`awk 'NR==1,NR==1 {print \$1}' MULTINODE`;
chomp($manynodes);
if($manynodes eq "Y") {
 $multinode=`awk 'NR==1,NR==1 {print \$0}' MULTINODEDATA`;
 chomp($multinode);
 $multinode=~ s/^\s+|\s+$//g;
}else{
return;
}

if ( -s $fin) {

open(IN,"<$fin");
open(TMP,">tmpfile");
$count=0;
while (<IN>) {
if (/bash/) {
}else{
$count=$count+1;
$line=$_;
chomp($line);
$line=~ s/^\s+|\s+$//g;
if($line eq "wait") {
print TMP "$line\n";
}else{
$replace=$count;
$thisnode="";
if ($multinode ne "") {
$thisnode=$multinode;
$thisnode =~ s/node/$replace/g;
}
$afin=join("",$fin,$replace);
$newline=join(" ",$thisnode,$afin,"&");
print TMP "$newline\n";
$line=~ s/&//;
if ("$fin" eq "runchg" || "$fin" eq "runpar" ) {
system("cat LOADQCP > $afin");
}
if ("$fin" eq "rundal") {
system("cat LOADDP > $afin");
}
open(F2,">>$afin");
print F2 "$line\n";
close F2;
system("chmod +x $afin");
}

}
}
close IN;
close TMP;
}
system("mv tmpfile $fin");

}



sub runallpar {
system("clear");
system("stty sane");


system("rm -f OUTLIST*");

 $filename='READY2GO';
 open($fh,'<:encoding(UTF-8)',$filename)
 or die "could not open '$filename' $!";
 while ($ready = <$fh>){
 chomp $ready;
 $ready2=&trim($ready);
 }
 close $fh;

 if ($ready2 < "2") {
 print "Before a job can be submitted (or re-submitted), SMFA must check the input\n";
 print "and perform the fragmentation.\n";
 print "SMFA will perform those preliminary steps now.\n";
 system("sleep 5");
 &fragonly;
 }
 system("echo 0 > READY2GO");

$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
$n++;
}
$nlast=$n-1;
open(TMP,"<IN_JOBTYPE");
while (<TMP>) {
if ($. == 2) {$deriv="$_"};
}
chomp $deriv;
$deriv=~ s/^\s+|\s+$//g;

$qmprog=$store[1];

# modify filenames in INLIST etc
if ( -s "INLISTCHG" ) {
open (IN,"<INLISTCHG");
open(TMP,">tempfile");
while (<IN>) {
$file="$_";
chomp($file);
$file=~ s/^\s+|\s+$//g;
@line=split(".grp",$file);
print TMP "$line[0].com\n";
}
close IN;
close TMP;
system("cp tempfile INLISTCHG");
system("cp INLISTCHG INLISTCHG_original");
}

#if ($qmprog == 1 || $qmprog == 3) {
#open (IN,"<INLIST");
#open(TMP,">tempfile");
#while (<IN>) {
#$file="$_";
#chomp($file);
#$file=~ s/^\s+|\s+$//g;
#@line=split(".com",$file);
#if ($qmprog == 1) {
#print TMP "$line[0].inp\n";
#}
#if ($qmprog == 3) {
#print TMP "$line[0].nw\n";
#}
#}
#}
# finished mods to filenames


&createrunscripts;

if( $deriv == 3 || $deriv == 4 ) {
&runoptpar;
exit(0);
}

if( $deriv == 5 ) {
&runscanpar;
exit(0);
}


if ($qmprog == 1) {
system("cp pbsfile_bare subrunallgam");
#system("cat LOADQCP >> subrunallgam");
system("cat LOADQCP > scr_gam");
system("cat $CODEDIR/scr_gam1 >>scr_gam");
system("cat LOADDP >>scr_gam");
system("cat $CODEDIR/scr_gam2 >>scr_gam");
system("cat ./scr_gam >> subrunallgam");
#system("cat $CODEDIR/scr_gam >> subrunallgam");
system("chmod +x subrunallgam");
system("qsub subrunallgam");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}

if ($qmprog == 2) {
system("cp pbsfile subrunallgau");
system("cat $CODEDIR/scr_gau >> subrunallgau");
system("chmod +x subrunallgau");
system("qsub subrunallgau");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}

if ($qmprog == 3) {
system("cp pbsfile_bare subrunallnwc");
#system("cat LOADQCP >> subrunallgam");
system("cat LOADQCP > scr_nwc");
system("cat $CODEDIR/scr_nwc1 >> scr_nwc");
system("cat LOADDP >> scr_nwc");
system("cat $CODEDIR/scr_nwc2 >> scr_nwc");
system("cat ./scr_nwc >> subrunallnwc");
#system("cat $CODEDIR/scr_nwc >> subrunallnwc");
system("chmod +x subrunallnwc");
system("qsub subrunallnwc");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}

if ($qmprog == 4) {
system("cp pbsfile subrunallqch");
system("cat $CODEDIR/scr_qch >> subrunallqch");
system("chmod +x subrunallqch");
system("qsub subrunallqch");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}


}

sub largest {

$qmprog=`awk 'NR==1,NR==1 {print \$0}' IN_PACKAGE`;
chomp($qmprog);

$frag=`awk '/for fragment number/ {print \$6}' OUT_ELECTRONS_SUMMARY`;
chomp($frag);

    $nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;
$deriv=`awk 'NR==2,NR==2 {print \$1}' IN_JOBTYPE`;
chomp $deriv;
$deriv=~ s/^\s+|\s+$//g;

if($qmprog == 1) {
$nchgs = 0;
if ($deriv == 0) {$prop = "RUNTYP=ENERGY"};
if ($deriv == 1) {$prop = "RUNTYP=GRADIENT"};
if ($deriv == 2) {$prop = "RUNTYP=HESSIAN"};

# get CONTRL group commands
 open(TMP,"<GAMcommands");
 $ic=0;
 while(<TMP>) {
 my $line="$_";
 chomp $line;
 $ic=$ic+1;
 $prop=join(" ",$prop,$line);
 }
  $ans=`awk 'NR==1,NR==1 {print \$1}' GAMbasis`;
  chomp $ans;
  my $lf=`wc -l GAMbasis`;
  chomp $lf;
  my @kf=split(/\s/,"$lf");
  $mf=$kf[0];
  $mf=$mf-1;
  system("tail -$mf GAMbasis > gbfile");
  $ch=`awk 'NR==1,NR==1{print \$1}' COORD$frag`;
  $mult=`awk 'NR==1,NR==1{print \$3}' COORD$frag`;
  chomp $ch;
  chomp $mult;
  open(TMP,">tmpfile");
  print TMP " \$CONTRL $prop MULT=$mult ICHARG=$ch \$END\n";
  close TMP;
  if ( -s "GAMextras") {system("cat GAMextras >> tmpfile")};
  open(TMP,">>tmpfile");
  print TMP " \$DATA\n";
  print TMP " GAMESS job\n";
  print TMP "c1\n";

  my $lf=`wc -l COORD$frag`;
  chomp $lf;
  my @kf=split(/\s/,"$lf");
  $mf=$kf[0];
  $mf=$mf-1;
  system("tail -$mf COORD$frag > tfile");
  system("$EXEDIR/convertcoords4games < tfile > gamesscoords");
  open(TMP1,"<gamesscoords");

  while (<TMP1>) {
   my $line1=$_;
   chomp $line1;
   my @bits=split(/\s+/,$line1);
   $ele=$bits[0];
   print TMP "$line1\n";
if ($ans eq "N") {
  open(GAMB,"<gbfile");
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
  system("cat gbfile >> tmpfile");
  open(TMP,">>tmpfile");
}
}
  print TMP "\$END\n";
  close TMP;
  close TMP1;
  system("mv tmpfile LARGE.inp");
}


if($qmprog == 2) {
system("cat IN_G09 > LARGE.com");
system("echo 'Bonded fragment' >> LARGE.com");
system("echo ' ' >> LARGE.com");
$file="COORD$frag";
system("cat $file >> LARGE.com");
system("echo ' ' >> LARGE.com");
}

if($qmprog == 3) {
  open(TMP,">tmpfile");
  my $my_coord = "COORD$frag";
  print TMP "START LARGE\n";
  print TMP "ECHO\n";
  $ch=`awk 'NR==1,NR==1{print \$1}' COORD$frag`;
  $mult=`awk 'NR==1,NR==1{print \$3}' COORD$frag`;
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
  my $lf=`wc -l COORD$frag`;
  chomp $lf;
  my @kf=split(/\s/,"$lf");
  $mf=$kf[0];
  $mf=$mf-1;
  close TMP;
  system("tail -$mf COORD$frag >> tmpfile");
  open(TMP,">>tmpfile");
  print TMP "END\n";
  print TMP "BASIS\n";
  close TMP;
  system("cat NWCbasis >> tmpfile");
  open(TMP,">>tmpfile");
  print TMP "END\n";
 open(FH,"<NWCcommands");
 while (<FH>) {
  if( $mult == 2 ) {
   s/SINGLET/$mform/;
  }
  print TMP "$_";
 }
  close TMP;
  system("mv tmpfile LARGE.nw");
}

if($qmprog == 4) {
system("cp IN_QCH LARGE.com");
system("echo '\$molecule' >> LARGE.com");
system("cat COORD$frag >> LARGE.com");
system("echo '\$end' >> LARGE.com");
}

system("echo ' ' >> OUT_SMFA");
system("echo 'An ab initio input file for this fragment is contained in' >> OUT_SMFA");
system("echo 'LARGE.com or LARGE.inp or LARGE.nw' >> OUT_SMFA");
system("echo ' ' >> OUT_SMFA");
}
        
sub runoptpar {


if($qmprog == 1) {
system("cp pbsfile_bare subrunoptgam");
system("cat LOADQCP >> subrunoptgam");
system("echo $CODEDIR/optscript_par >> subrunoptgam");
system("chmod +x subrunoptgam");
system("qsub subrunoptgam");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}

if($qmprog == 2) {
system("cp pbsfile subrunoptgau");
system("echo $CODEDIR/optscript_par >> subrunoptgau");
system("chmod +x subrunoptgau");
system("qsub subrunoptgau");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}

if($qmprog == 3) {
system("cp pbsfile_bare subrunoptnwc");
system("cat LOADQCP >> subrunoptnwc");
system("echo $CODEDIR/optscript_par >> subrunoptnwc");
system("chmod +x subrunoptnwc");
system("qsub subrunoptnwc");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}

if($qmprog == 4) {
system("cp pbsfile subrunoptqch");
system("echo $CODEDIR/optscript_par >> subrunoptqch");
system("chmod +x subrunoptqch");
system("qsub subrunoptqch");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}

}

sub runscanpar {

if($qmprog == 1) {
system("cp pbsfile_bare subrunscangam");
system("cat LOADQCP >> subrunscangam");
system("echo $CODEDIR/scanscript_par >> subrunscangam");
system("chmod +x subrunscangam");
system("qsub subrunscangam");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}

if($qmprog == 2) {
system("cp pbsfile subrunscangau");
system("echo $CODEDIR/scanscript_par >> subrunscangau");
system("chmod +x subrunscangau");
system("qsub subrunscangau");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}

if($qmprog == 3) {
system("cp pbsfile_bare subrunscannwc");
system("cat LOADQCP >> subrunscannwc");
system("echo $CODEDIR/scanscript_par >> subrunscannwc");
system("chmod +x subrunscannwc");
system("qsub subrunscannwc");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}

if($qmprog == 4) {
system("cp pbsfile subrunscanqch");
system("echo $CODEDIR/scanscript_par >> subrunscanqch");
system("chmod +x subrunscanqch");
system("qsub subrunscanqch");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);
}





}


sub restartpar {

open(SM,">>OUT_SMFA");
print SM "\n";
print SM "Incomplete calculations have been restarted.\n";
print SM "So, YOU SHOULD IGNORE THE RESULT ABOVE.\n";
print SM "The completed result will be added below.\n";
print SM "\n";
close SM;

if ( -s "INLIST_original" ) {
 system("cp INLIST_original INLIST");
}
if ( -s "INLISTDAL_original" ) {
 system("cp INLISTDAL_original INLISTDAL");
}

if (($qmprog == 1) || ($qmprog == 3)) {
 system("cp pbsfile_bare subrestartallpar");
 system("cat LOADQCP >> subrestartallpar");
}else{
system("cp pbsfile subrestartallpar");
}
system("cat $CODEDIR/restartpar >> subrestartallpar");
system("chmod +x subrestartallpar");
system("qsub subrestartallpar");
print "\n";
print "Ab initio jobs submitted\n";
print "\n";
print "The results will be appended to OUT_SMFA\n";
#system("sleep 4");
system("stty sane");
#system("clear");
exit(0);


sub makeqmdal {

open(FI,">pbsfile_bare1");
open(PB,"<pbsfile");
open(QCP,">LOADQCP");

 my $skip=0;
 my $writit=0;
 while(<PB>) {
  if(/that make the quantum chemistry package/) {
   $skip=$. + 2;
  }
  if ($. == $skip) {$writit=1};
  if ($_ eq "\n") {$writit=0};
  if( $writit == 1) {
   print QCP "$_";
  }else{
   print FI "$_";
  }
  }

close PB;
close QCP;
close FI;

open(FI,">pbsfile_bare");
open(PB,"<pbsfile_bare1");
open(DP,">LOADDP");

 my $skip=0;
 my $writit=0;
 while(<PB>) {
  if(/that make the DALTON  package/) {
   $skip=$. + 2;
  }
  if ($. == $skip) {$writit=1};
  if ($_ eq "\n") {$writit=0};
  if( $writit == 1) {
   print DP "$_";
  }else{
   print FI "$_";
  }
 }
close PB;
close DP;
close FI;
system("rm pbsfile_bare1");

}
}

sub trimxyz {
$targetfile=$_[0];

open(IN,"<$targetfile");
open(OUT,">tmpfile");
while (<IN>) {
if("$_" ne "\n") {
$line=&ltrim($_);
print OUT "$line";
}else{
print OUT "\n";
}
}
close IN;
close OUT;
system("cp tmpfile $targetfile");
system("rm tmpfile");
}

sub getncpus {
$line=`awk /NCPUS=/ pbsfile`;
if($line eq "\n" || $line eq "") {
 $line=`awk /ncpus=/ pbsfile`;
}
if($line eq "\n" || $line eq "") {
 print "The pbsfile is missing a statement containing 'ncpus='\n";
 print "SMFA will terminate\n";
 system("sleep 5");
 exit(0);
}
}


