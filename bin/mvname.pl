#!/usr/bin/perl

# usage is mvname INLIST    or mvname INLISTCHG

$list=$ARGV[0];
chomp($list);
$list=~ s/^\s+|\s+$//g;

$package =`awk 'NR==1,NR==1 {print \$0}' IN_PACKAGE`;
chomp($package);

if( $list eq "INLIST") {
 $line=`awk /ncpus=/ pbsfile`;
 @bline=split("=",$line);
 $ncpus=$bline[1];
 chomp($ncpus);
 $ncpus=~ s/^\s+|\s+$//g;
}

if ( $list eq "INLISTCHG") {
 $ncpus=`awk 'NR==1,NR==1 {print \$1}' NCPUSCH`;
 chomp($ncpus);
 $ncpus=~ s/^\s+|\s+$//g;
}

for ($i=1;$i <= $ncpus;$i++) {
$outnum=join("","OUT",$list,$i);
open(OUTLOC,"<$outnum");
while (<OUTLOC>) {
$job=$_;
chomp($job);
$job=~ s/^\s+|\s+$//g;
$del=join("","_",$i);
@st=split($del,$job);
$stem=$st[0];
$com=join("",$stem,".com");
$comloc=join("",$stem,$del,".com");
$log=join("",$stem,".log");
system("mv $job $log");

if ( -s $comloc) {
system("mv $comloc $com");
}

if ($package == 1) {
$inp=join("",$stem,".inp");
$inploc=join("",$stem,$del,".inp");
if ( -s $inploc) {
system("mv $inploc $inp");
}
$dat=join("",$stem,$del,".dat");
$dat0=join("",$stem,".dat");
if (-s $dat) {
system("mv $dat $dat0");
}
}

if ($package == 2) {
 $fchkloc=join("",$stem,$del,".fchk");
 $fchk=join("",$stem,".fchk");
 if( -s "$fchkloc") {
  system("mv $fchkloc $fchk");
 }
}

if ($package == 3) {
 $nw=join("",$stem,$del,".nw");
 $nw0=join("",$stem,".nw");
 system("mv $nw $nw0");
 $esp=join("",$stem,$del,".esp");
 $esp0=join("",$stem,".esp");
 if ( -s "$esp") {
  system("mv $esp $esp0");
 }
 $hes=join("",$stem,$del,".hess");
 $hes0=join("",$stem,".hess");
 if ( -s "$hes") {
  system("mv $hes $hes0");
 }
}

if ($package == 4) {
 $fchkloc=join("",$stem,$del,".com",".fchk");
 $fchk=join("",$stem,".com",".fchk");
 if( -s "$fchkloc") {
  system("mv $fchkloc $fchk");
 }
}
# end the while
}
close OUTLOC;
#end the for i loop
}


