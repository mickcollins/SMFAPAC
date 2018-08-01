#!/usr/bin/perl

$snum=$ARGV[0];
chomp($snum);
$snum=~ s/^\s+|\s+$//g;

&getcombdata("$snum");

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



