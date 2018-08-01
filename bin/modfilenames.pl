#!/usr/bin/perl

&modfilenames;


sub modfilenames {

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
}




