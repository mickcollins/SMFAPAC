#!/usr/bin/perl -w

$check="";
$new="";
while ( @ARGV ) {
  $flag=shift;
  $check=shift, next if $flag eq "-check";
  $new=shift, next if $flag eq "-new";
  $check=$flag, next if $check eq "";
  $new=$flag, next if $new eq "";
  die "Don't understand command-line item $flag. Stopped";
}

open (CHECK,"$check") or die "Can't open file $check. Stopped";
open (NEW,"$new") or die "Can't open file $new. Stopped";

$mismatches=0;
$newlines="";
$checklines="";
while ( <NEW> ) {
  s/\r//;
  next if /^ *$/;
  while ( $c=<CHECK> ) {
    last unless $c=~/^ *$/;
  }
  next if /^Starting at/ || /^Finished at/ || /^CPU time used/;
  s/-(0\.0+)\b/ $1/g;
  $c=~s/-(0\.0+)\b/ $1/g;
  if ($_ ne $c) {
    # print $_, $c;
    $newlines.=$_;
    $checklines.=$c;
    print "New output file $new differs from check file $check.\n" unless $mismatches;
    $mismatches++;
  }
  elsif ( $newlines ) {
    print "New file:\n$newlines\nCheck file:\n$checklines\n\n";
    $newlines="";
    $checklines="";
  }
  if ( $mismatches>20 ) {
    print "Too many differences\n";
    exit 1;
  }
}

$c=<CHECK>;
print "New output file ended before check file" unless eof(CHECK);

close NEW;
close CHECK;
print "Checked -- O.K.\n" unless $mismatches;

