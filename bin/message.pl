#!/usr/bin/perl

# this script looks for error messages about "failed" in OUT_SMFA
# and writes warnings/messages to OUT_SMFA

$errors=0;
$nchg=0;
open(IN,"<OUT_SMFA");
while (<IN>) {
if (/may have failed/) {
 $errors = $errors + 1;
 if(/charge/) {$nchg = $nchg + 1};
}
if (/SHOULD IGNORE THE RESULT/) {
 $errors = 0;
 $nchg=0;
}
}
close IN;

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
&rmfiles;
}

sub rmfiles {

$junk=`rm -f list.* 2> /dev/null`;
$junk=`rm -f tmpf* 2> /dev/null`;
$junk=`rm -f tfil* 2> /dev/null`;
$junk=`rm -f IN_G09_temp* 2> /dev/null`;
$junk=`rm -f IN_G09_current* 2> /dev/null`;
$junk=`rm -f polarfiles* 2> /dev/null`;
$junk=`rm Lev0* 2> /dev/null`;
$junk=`rm Lev1* 2> /dev/null`;
$junk=`rm *COORD* 2> /dev/null`;
$junk=`rm chL* 2> /dev/null`;
$junk=`rm chFR* 2> /dev/null`;
$junk=`rm chnb* 2> /dev/null`;
$junk=`rm chab* 2> /dev/null`;
$junk=`rm charge* 2> /dev/null`;
$junk=`rm OUT_NCOORD_FRAG_Lev0 2> /dev/null`;
$junk=`rm OUT_L1L1_data 2> /dev/null`;
$junk=`rm CHARGECOORDS 2> /dev/null`;
$junk=`rm OUT_CHARGEDGROUPS 2> /dev/null`;
$junk=`rm VIEWCHARGEDGROUPS 2> /dev/null`;
$junk=`rm OUT_NCOORD_FRAG 2> /dev/null`;
$junk=`rm fragnu* 2> /dev/null`;
$junk=`rm archive* 2> /dev/null`;
$junk=`rm thisenergy* 2> /dev/null`;
$junk=`rm theseforces* 2> /dev/null`;
$junk=`rm gamesscoords* 2> /dev/null`;
$junk=`rm thesecoords* 2> /dev/null`;
$junk=`rm thesedipdr* 2> /dev/null`;
$junk=`rm NearInt_energy* 2> /dev/null`;
$junk=`rm gbfile* 2> /dev/null`;
$junk=`rm gamjunk* 2> /dev/null`;
$junk=`rm abhessians* 2> /dev/null`;

$junk=`rm FRAG*com 2> /dev/null`;
$junk=`rm FRAG*movecs 2> /dev/null`;
$junk=`rm FRAG*db 2> /dev/null`;
$junk=`rm FRAG*nw 2> /dev/null`;
$junk=`rm FRAG*inp 2> /dev/null`;

$junk=`rm ab*com 2> /dev/null`;
$junk=`rm ab*nw 2> /dev/null`;
$junk=`rm ab*inp 2> /dev/null`;
$junk=`rm ab*movecs 2> /dev/null`;
$junk=`rm ab*db 2> /dev/null`;

$junk=`rm nb*com 2> /dev/null`;
$junk=`rm nb*inp 2> /dev/null`;
$junk=`rm nb*cart 2> /dev/null`;
$junk=`rm nb*esp 2> /dev/null`;
$junk=`rm nb*grid 2> /dev/null`;
$junk=`rm nb*mol 2> /dev/null`;
$junk=`rm nb*dal 2> /dev/null`;
$junk=`rm nb*.nw 2> /dev/null`;
$junk=`rm nb*.db 2> /dev/null`;
$junk=`rm nb*.mol 2> /dev/null`;
$junk=`rm nb*.movecs 2> /dev/null`;
$junk=`rm nb*.plt 2> /dev/null`;
$junk=`rm nb*.xyz 2> /dev/null`;
$junk=`rm nb*.q 2> /dev/null`;
$junk=`rm nb*polar* 2> /dev/null`;
$junk=`rm nb*Stone 2> /dev/null`;

$junk=`rm nb*punch 2> /dev/null`;
$junk=`rm nb*gdma* 2> /dev/null`;
$junk=`rm rub* 2> /dev/null`;
$junk=`rm grp* 2> /dev/null`;
$junk=`rm OUTINLIST* 2> /dev/null`;
$junk=`rm OUTINLISTCHG* 2> /dev/null`;
$junk=`rm fort* 2> /dev/null`;
$junk=`rm abfiles* 2> /dev/null`;
$junk=`rm ABbgidentities* 2> /dev/null`;
$junk=`rm OUT_ABGROUPS 2> /dev/null`;
$junk=`rm junk* 2> /dev/null`;
$junk=`rm OUT_ALLOC_FRAG 2> /dev/null`;
$junk=`rm OUT_CHARGE_CHARGE 2> /dev/null`;
$junk=`rm OUT_ELECTRONS 2> /dev/null`;
$junk=`rm OUT_ELECTRONS_SUMMARY 2> /dev/null`;
$junk=`rm OUT_FINAL_L1_DATA 2> /dev/null`;
$junk=`rm OUT_FRAGMENTS_Lev0 2> /dev/null`;
$junk=`rm OUT_FRAGMENTS_Lev1 2> /dev/null`;
$junk=`rm OUT_FRAGMENTS_LevX 2> /dev/null`;
$junk=`rm OUT_GROUPCONNECTIVITY 2> /dev/null`;
$junk=`rm OUT_L1L1_AB_data 2> /dev/null`;
$junk=`rm OUT_L1_FINALSIGNS 2> /dev/null`;
$junk=`rm OUT_LIST_FRAG 2> /dev/null`;
$junk=`rm OUT_METALGROUPS 2> /dev/null`;
$junk=`rm OUT_NGROUPS 2> /dev/null`;
$junk=`rm OUT_NONMETALCHARGES 2> /dev/null`;
$junk=`rm RAWDATA 2> /dev/null`;
$junk=`rm RUNDATA 2> /dev/null`;

$junk=`rm combinedABNBderivs 2> /dev/null`;
$junk=`rm combinedABbgderivs 2> /dev/null`;
$junk=`rm combinedFRAGbgderivs 2> /dev/null`;
$junk=`rm combinedFRAGderivs 2> /dev/null`;
$junk=`rm Dispderivatives 2> /dev/null`;
$junk=`rm Doublecount_derivs 2> /dev/null`;
$junk=`rm Electrostatic_derivs 2> /dev/null`;

$junk=`rm frags.out_Lev0 2> /dev/null`;
$junk=`rm frags.out_Lev1 2> /dev/null`;
$junk=`rm nbfilesN* 2> /dev/null`;
$junk=`rm signs.out_Lev0 2> /dev/null`;
$junk=`rm signs.out_Lev1 2> /dev/null`;
$junk=`rm INLIST.* 2> /dev/null`;
$junk=`rm INLISTDAL.* 2> /dev/null`;

$junk=`rm qchinput* 2> /dev/null`;
$junk=`rm getdata* 2> /dev/null`;

}



