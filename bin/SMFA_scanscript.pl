#!/usr/bin/perl

#use warnings;
use Config;
use lib "CPAN";
use Shell qw(date pwd);
use File::Copy 'cp';
use Cwd;

$CODEDIR=`awk 'NR==1,NR==1 {print \$1}' CODENAME`;
chomp($CODEDIR);
$CODEDIR=~ s/^\s+|\s+$//g;

$EXEDIR=`awk 'NR==1,NR==1 {print \$1}' EXENAME`;
chomp($EXEDIR);
$EXEDIR=~ s/^\s+|\s+$//g;

&scan;
&getscanenergies;

exit(0);

