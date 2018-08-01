#!/usr/bin/perl

$CODEDIR=`awk 'NR==1,NR==1 {print \$1}' CODENAME`;
chomp($CODEDIR);
$CODEDIR=~ s/^\s+|\s+$//g;

$EXEDIR=`awk 'NR==1,NR==1 {print \$1}' EXENAME`;
chomp($EXEDIR);
$EXEDIR=~ s/^\s+|\s+$//g;

&createrunscripts;

sub createrunscripts {

$line=`awk /ncpus=/ pbsfile`;
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
&multnodeinput("gauchinput");
}
if ($package == 3) {
&multnodeinput("nwcinput");
&multnodeinput("nwcchinput");
}
if ($package == 4) {
&multnodeinput("qchinput");
&multnodeinput("qchchinput");
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
$d="d";
$del="1,$ncpus$d";
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
$del=join("","1,",$ncpuschdal,"d");
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


