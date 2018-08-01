#!/usr/bin/perl

$CODEDIR=`awk 'NR==1,NR==1 {print \$1}' CODENAME`;
chomp($CODEDIR);
$CODEDIR=~ s/^\s+|\s+$//g;

$EXEDIR=`awk 'NR==1,NR==1 {print \$1}' EXENAME`;
chomp($EXEDIR);
$EXEDIR=~ s/^\s+|\s+$//g;

&dodaltondisp;




sub extract_disp_dal  {
 system("rm PAA.out");
 open($PAA_output,">PAA.out");

 my $ngroups = `ls grp.*.mol | wc -l`;
  for(my $group=1;$group<=$ngroups;$group++)
  {
   my $filebase="grp.$group";
   $disp_out="$filebase.disp";
   $dalton_input="$filebase.mol"; 

      open($PAA_In, ">$disp_out");

      open($DAL_IN,"$dalton_input") or die "Unable to open $DAL_IN";
      @dalton_in = <$DAL_IN>;
      while (my $DLine = shift (@dalton_in))
      {
         if ($DLine =~ /^(\w+\s+)((-)*\d+\.\d+\s+)+/) #match  line that looks like coordinates
         {
            $num_atoms+=1;
            $mol_coords[$num_atoms]=$DLine;
         }
      }
      close $DAL_IN ;
      print $PAA_In "$num_atoms\n";
     for (my $ii=1; $ii<=$num_atoms; $ii++)
     {
      print $PAA_In $mol_coords[$ii];
     }
      undef(@mol_coords);

      @dalton_output=`grep -A 17 -B 1 "GRIDSQ" $filebase.out`;

      foreach $line (@dalton_output)
      {
         if ($line !~ /(GRIDSQ|--)/) # the "--" is output by grep wen there's a blank line
              {
                      print $PAA_In "$line";
              }
      }
      #Now get the PAA value and print it to a file
      @PAA_guff=`$EXEDIR/PAA <$filebase.disp`;
      my $last=$num_atoms+3;
      print $PAA_output "$group", "   ",  "$PAA_guff[$last]";
      $num_atoms=0;
  }
}

sub ReadPolar_dal {
system("rm Polarisabilities.out");
open(my $OUT, ">Polarisabilities.out");

my $nfrgs = `ls nb.*.0-polar.out | wc -l`;
chomp $nfrgs;
for (my $ic=1; $ic<=$nfrgs; $ic++) {
$job="nb.$ic.0-polar.out";
&getpolar_dal("$job");
$iso=($Polar[1] + $Polar[5] + $Polar[9])/3;
print $OUT "$ic\t$Polar[1]\t$Polar[4]\t$Polar[5]\t$Polar[7]\t$Polar[8]\t$Polar[9]\t$iso\n";

}
close $OUT;
}

sub getpolar_dal  {
 $file="<$_[0]";
 $ans=`grep 'FREQUENCY INDEPENDENT SECOND ORDER PROPERTIES' $file | wc -l`;
 chomp($ans);
 if ( $ans == 1) {
  $Polar[1]=`awk '/XDIPLEN  ; XDIPLEN/ {print \$8}' $file`;
  $Polar[2]=`awk '/XDIPLEN  ; YDIPLEN/ {print \$8}' $file`;
  $Polar[3]=`awk '/XDIPLEN  ; ZDIPLEN/ {print \$8}' $file`;
  $Polar[4]=$Polar[2];
  $Polar[5]=`awk '/YDIPLEN  ; YDIPLEN/ {print \$8}' $file`;
  $Polar[6]=`awk '/YDIPLEN  ; ZDIPLEN/ {print \$8}' $file`;
  $Polar[7]=$Polar[3];
  $Polar[8]=$Polar[6];
  $Polar[9]=`awk '/ZDIPLEN  ; ZDIPLEN/ {print \$8}' $file`;
  for ($i=1;$i <= 9;$i++) {
   chomp($Polar[$i]);
  }
  return;
 }

 open(FRAGJOB,"<$_[0]");
 $skip=0;
 $skip2=0;
 $writit=0;
 $nvar=1;
while (<FRAGJOB>)  {
 if(/Static polarizabilities \(au\)/) {
   $skip=$. + 5;
   $skip2=$. + 8;
}
if ($. == $skip) {$writit=1};
if ($. == $skip2) {$writit=0};
      if( $writit == 1) {
       $line1="$_";
       chomp $line1;
       @line=split(/\s+/,"$line1");
       $Polar[$nvar]=$line[2];
       $nvar=$nvar+1;
       $Polar[$nvar]=$line[3];
       $nvar=$nvar+1;
       $Polar[$nvar]=$line[4];
       $nvar=$nvar+1;
      }
 }
 close FRAGJOB;

}


sub dodaltonpolar {
$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
my $n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
$n++;
}
$n=$n-1;
close $fh;

# find which package
# and what electronic structure method
if($store[1] == 3) {
#its nwc
$here=`grep "TASK SCF" "ABSETTINGS" | wc -l`;
chomp $here;
if($here > 0) {
 $method="HF";
}else{
 $method="MP2";
}
}

if($store[1] == 1) {
#its gam
$here=`grep "HF" "ABSETTINGS" | wc -l`;
chomp $here;
if($here > 0) {
 $method="HF";
} 
$here=`grep "MPLEVL" "ABSETTINGS" | wc -l`;
chomp $here;
if($here > 0) {
 $method="MP2";
}
}


for (my $ic=0; $ic<=$n; $ic++) {
 if ("$store[$ic]" eq "Is long range dispersion accounted for?")  {
  $ie=$ic+1;
  $disp=$store[$ie];
 }
}

# make the input files that DALTON needs
&daltonpolar;

$nchgs =`awk 'NR==2,NR==2 {print \$1}' IN_CHARGES`;
chomp $nchgs;
if ($store[1] == 1) {$nchgs=0};

if ("$disp" eq "Y") {
&SMFAmkdalpolar;
#&SMFArundalpolar;
 } else {
   if ($nchgs eq 0) {
&SMFAmkdalpolar;
#&SMFArundalpolar;
   }
}

}


sub SMFAmkdalpolar {
system("ls -1 nb.*.0-polar.com > polarfiles");
$filename='polarfiles';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";

open(IN,">>INLISTDAL");

while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 system("$EXEDIR/convert2dalton < nb.$ic1.0-polar.com > nb.$ic1.0-polar.mol");
 print IN "nb.$ic1.0-polar.mol\n";
 $mult=`awk 'NR==1,NR==1 {print \$2}' nb.$ic1.0-polar.com`;
 chomp $mult;
 if($method eq "HF") {
  if($mult == 1) {
   system("cp IN_DAL_POLAR_HF_S nb.$ic1.0-polar.dal");
   }else{
   system("cp IN_DAL_POLAR_HF_D nb.$ic1.0-polar.dal");
   }
  }else{
    if($mult == 1) {
   system("cp IN_DAL_POLAR_MP2_S nb.$ic1.0-polar.dal");
   }else{
#   system("cp IN_DAL_POLAR_MP2_D nb.$ic1.0-polar.dal");
#   # approximation to stop dalton crashing with MP2 doublets
   system("cp IN_DAL_POLAR_HF_D nb.$ic1.0-polar.dal");
   }
  }
}
close $fh;
close IN;
}


sub SMFArundalpolar {

 system("ls -1 nb*.0-polar.mol > nbfilesp");
$filename='nbfilesp';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
while ($row = <$fh>){
 chomp $row;
 my @fields=split(/\./,$row);
 my $ic1=$fields[1];
 system("dalton -N 1 -o nb.$ic1.0-polar.out nb.$ic1.0-polar");
 }
close $fh;
}

sub daltonpolar {

open(TMP,">IN_DAL_POLAR_HF_D");
print TMP "**DALTON INPUT\n";
print TMP ".RUN PROPERTIES\n";
print TMP "**WAVE FUNCTIONS\n";
print TMP ".HF\n";
print TMP "*CONFIGURATION INPUT\n";
print TMP ".SPIN MULTIPLICITY\n";
print TMP "  2\n";
print TMP "**PROPERTIES\n";
print TMP ".POLARI\n";
print TMP "**END OF DALTON INPUT\n";
close TMP;

open(TMP,">IN_DAL_POLAR_HF_S");
print TMP "**DALTON INPUT\n";
print TMP ".RUN PROPERTIES\n";
print TMP "**WAVE FUNCTIONS\n";
print TMP ".HF\n";
print TMP "**PROPERTIES\n";
print TMP ".POLARI\n";
print TMP "**END OF DALTON INPUT\n";
close TMP;

open(TMP,">IN_DAL_POLAR_MP2_D");
print TMP "**DALTON INPUT\n";
print TMP ".RUN RESPONSE\n";
print TMP "**WAVE FUNCTIONS\n";
print TMP ".HF\n";
print TMP ".MP2\n";
print TMP "*CONFIGURATION INPUT\n";
print TMP ".SPIN MULTIPLICITY\n";
print TMP "  2\n";
print TMP "**RESPONSE\n";
print TMP ".SOPPA\n";
print TMP "*LINEAR\n";
print TMP ".DIPLEN\n";
print TMP ".FREQUENCIES\n";
print TMP " 1\n";
print TMP " 0.0\n";
print TMP "**END OF DALTON INPUT\n";
close TMP;

open(TMP,">IN_DAL_POLAR_MP2_S");
print TMP "**DALTON INPUT\n";
print TMP ".RUN RESPONSE\n";
print TMP "**WAVE FUNCTIONS\n";
print TMP ".HF\n";
print TMP ".MP2\n";
print TMP "**RESPONSE\n";
print TMP ".SOPPA\n";
print TMP "*LINEAR\n";
print TMP ".DIPLEN\n";
print TMP ".FREQUENCIES\n";
print TMP " 1\n";
print TMP " 0.0\n";
print TMP "**END OF DALTON INPUT\n";
close TMP;
}


sub dodaltondisp {
# make IN_DALTON
open(TMP,">IN_DALTON");
print TMP "basis=6-311G**\n";
print TMP "ATOMBASIS\n";
print TMP "molecule\n";
print TMP "Disp calculation, no charge field\n";
close TMP;
# make IN_DISP
open(TMP,">IN_DISP");
print TMP "**DALTON INPUT\n";
print TMP ".RUN RESPONSE\n";
print TMP "**WAVE FUNCTIONS\n";
print TMP ".HF\n";
print TMP "**RESPONSE\n";
print TMP ".MAXRM\n";
print TMP "400\n";
print TMP "*C6\n";
print TMP ".MAXMOM\n";
print TMP "20\n";
print TMP ".GSLEGN\n";
print TMP ".DIPLEN\n";
print TMP ".PRINT\n";
print TMP "10\n";
print TMP "**END OF DALTON INPUT\n";
close TMP;

$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
$n++;
}
close $fh;
$n=$n-1;

$qmprog=$store[1];
for (my $ic=0; $ic<=$n; $ic++) {
 if ("$store[$ic]" eq "Is long range dispersion accounted for?")  {
  $ie=$ic+1;
  $disp=$store[$ie];
 }
}

if("$disp" eq "Y") {

&SMFAdaldispinputs_MAC;
#&SMFArundispdal;
}

}


sub SMFAdaldispinputs_MAC  {

#open(IN,">>INLISTDAL");
 my $nfrgs = `ls grp.*.mol | wc -l`;
  for (my $ic=1; $ic<=$nfrgs; $ic++)
  {
   my $filebase="grp.$ic";
   &Daltonify("$filebase.mol");
   system ("ln -s IN_DISP grp.$ic.dal");
#   print IN "$filebase.mol\n";
  }
close IN;
}

sub Daltonify  {

 $filename="$_[0]";

$pwd=`pwd`;
chomp $pwd;

$header="IN_DALTON";

open (my $comfile, "$filename") or die "Unable to open $filename";
@coords=<$comfile>;
close $comfile;

open (my $comfile, ">$filename") or die "Unable to open $filename for writing";
read_dump_headers($header,$comfile);
if (defined $basis)
{
   dump_coords($comfile,$basis,$H_basis);
}
else
{
   print $comfile @coords;
}
   print $comfile "\n";

close $comfile;
}

sub read_dump_headers 
{
        $Headerfile = $_[0];
             $fh = $_[1];

         if ($filename =~ /mol$/){    #dalton molecule files always end in .mol
             $filename =~ s/\.[^.]*$//;
         }

        open(HEADER,"$Headerfile") or die "Unable to open $Headerfile";
        @HeaderLines = <HEADER>;
        while ($HLine = shift (@HeaderLines))
        {
                if ($HLine =~ /molecule/) #find the title line
                {
                        print $fh "$filename\n";
                }
                elsif($HLine =~ /^H(B|b)asis=/) # stick a line with the basis in this file
                {
                   $H_basis=$';
                   chomp $H_basis;
                }
                elsif($HLine =~ /^(B|b)asis=/) # stick a line with the basis in this file
                {
                   $basis=$';
                   chomp $basis;
                }
                else
                {
                        print $fh "$HLine";
                }
        }
        close HEADER;
        if(!defined($H_basis))
        {
           #If we haven't defined a separate basis set for hydrogen atoms, then hydrogen atoms use the same
           #           #basis as all other atoms
           $H_basis=$basis
        }
}

sub dump_coords
{
             $fh = $_[0];
        $basis=$_[1];
        $H_basis=$_[2];

        while ($HLine = shift (@coords))
        {
           #print $HLine;
           if ($HLine =~ /Basis/) #find the Basis line 
              {
                 $tmp_start_of_line=$`;
                 if($tmp_start_of_line =~ /^Charge=1.0/) #Line refers to hydrogen
                 {
                    print $fh "$tmp_start_of_line Basis=$H_basis\n";
                 }
                 else
                 {
                   print $fh "$tmp_start_of_line Basis=$basis\n";
                }
              }
              else
              {
                 #           print $fh;
                      print $fh "$HLine";
              }
        }

}

sub SMFArundispdal  {

system("rm -f *.dal");

 my $nfrgs = `ls grp.*.mol | wc -l`;
  for (my $ic=1; $ic<=$nfrgs; $ic++)
  {
   system ("ln -s IN_DISP grp.$ic.dal");
   system("dalton -N 1 -o grp.$ic.out grp.$ic");
  }
}

