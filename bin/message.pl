#!/usr/bin/perl

# this script looks for error messages about "failed" in OUT_SMFA
# and writes warnings/messages to OUT_SMFA

$errors=`grep 'may have failed' OUT_SMFA | wc -l`;
$errors=~ s/^\s+|\s+$//g;
system("grep 'may have failed' OUT_SMFA > junk1");


$nchg=`grep charge junk1 | wc -l`;
$nchg=~ s/^\s+|\s+$//g;

if($nchg > 0) {
system("grep charge junk1 > junk2");
open(SM,">>OUT_SMFA");
print SM "IT APPEARS THAT A 'charge' CALCULATION HAS FAILED.\n";
print SM "To verify this, you should look at the output 'log' files,\n";
print SM "associated with the corresponding input files:\n";
close SM;
system("cat junk2 >> OUT_SMFA");
open(SM,">>OUT_SMFA");
print SM "If 'charge' jobs have failed, YOU CANNOT USE RESTART,\n";
print SM "but must simply run the whole job as usual\n";
print SM " (having corrected any problem created by the user).\n";
close SM;
}

if ( $errors != 0 && $nchg == 0 ) {
open(SM,">>OUT_SMFA");
print SM "IT APPEARS THAT SOME AB INITIO CALCULATION HAS FAILED.\n";
print SM "To verify this, you should look at the output 'log' files,\n";
print SM "associated with the corresponding input files:\n";
close SM;
system("cat junk1 >> OUT_SMFA");
open(SM,">>OUT_SMFA");
print SM "If all jobs have completed correctly, then the output above is correct.\n";
print SM "If some jobs have indeed failed, BUT THE INPUT IS CORRECT, then failure\n";
print SM "may have been due to a system/hardware problem.\n";
print SM "IN THIS CASE, YOU CAN SIMPLY USE THE RESTART FACILITY FOR PARALLEL JOBS\n";
close SM;
}

if( $errors == 0 && $nchg == 0 ) {
system("rm -f *.chk");
system("rm -f *.fchk");
#&rmfiles;
}

sub rmfiles {

system("rm -f list.*");
system("rm -f tmpf*");
system("rm -f tfil*");
system("rm -f IN_G09_temp*");
system("rm -f IN_G09_current*");
system("rm -f polarfiles*");
system("rm Frag*");
system("rm Lev0*");
system("rm Lev1*");
system("rm *COORD*");
system("rm chL*");
system("rm chFR*");
system("rm chnb*");
system("rm chab*");
system("rm charge.*");
system("rm OUT_NCOORD_FRAG_Lev0");
system("rm OUT_L1L1_data");
system("rm CHARGECOORDS");
system("rm OUT_CHARGEDGROUPS");
system("rm VIEWCHARGEDGROUPS");
system("rm OUT_NCOORD_FRAG");
system("rm fragnu*");

}

