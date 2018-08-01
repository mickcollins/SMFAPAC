#!/usr/bin/perl
# this is the runlist program

$job=$ARGV[0];
chomp($job);
$job=~ s/^\s+|\s+$//g;
$list=$ARGV[1];
chomp($list);
$list=~ s/^\s+|\s+$//g;

$num=$ARGV[2];
chomp($num);
$num=~ s/^\s+|\s+$//g;

$sl=$num/100;

open(RUN,">>RUNNING");
print RUN "$job\n";

@line=split("IN",$list);
$fin=$line[1];
$out="OUT$fin";

system("rmdir flagdr 2> /dev/null");

open(DON,">>$out");

$nproc0=1;

@line=split(".mol",$job);
$stem=$line[0];
$pol=index($stem,"polar");
$old=`grep $stem.out $out | wc -l`;
chomp($old);
$old=~ s/^\s+|\s+$//g;
if ( $old == 0 ) {
$junk=`dalton -N $nproc0 -o $stem.out $stem`;
#system("dalton -N $nproc0 -o $stem.out $stem");
 $ok1=`grep 'End of Dynamic Property Section' $stem.out | wc -l`;
 $ok2=`grep 'End of Static Property Section' $stem.out | wc -l`;
 chomp($ok1);
 chomp($ok2);
 $ok=$ok1 + $ok2;
 if ( $ok == 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation may have failed for $stem.mol\n";
  close TMP;
 }else{
  print DON "$stem.out\n";
 }
}


$go=0;
$end=1;
while ($go < $end ) {
RETRY:
system("mkdir flagdr 2> /dev/null");
if( $? > 8 ) {
system("sleep 0.028");
goto RETRY;
}else{
system("sed -i.$num '1d' $list");
system("rmdir flagdr 2> /dev/null");
}
$job =`awk 'NR==1,NR==1 {print \$1}' $list.$num`;
chomp($job);
$job=~ s/^\s+|\s+$//g;
 if ( "$job" eq "" ) {
  $go=2;
  last;
 }

 @line=split(".mol",$job);
 $stem=$line[0];
 chomp($stem);
 $stem=~ s/^\s+|\s+$//g;
 $old=`grep $stem.out $out | wc -l`;
 chomp($old);
 $old=~ s/^\s+|\s+$//g;
#system("sleep 0.028");
 $thisone=`grep $job RUNNING | wc -l`;
 chomp($thisone);
 $thisone=~ s/^\s+|\s+$//g;
 if ( $thisone == 0 && $old == 0 ) {
 print RUN "$job\n";
$junk=`dalton -N $nproc0 -o $stem.out $stem`;
#system("dalton -N $nproc0 -o $stem.out $stem");
 $ok1=`grep 'End of Dynamic Property Section' $stem.out | wc -l`;
 $ok2=`grep 'End of Static Property Section' $stem.out | wc -l`;
 chomp($ok1);
 chomp($ok2);
 $ok=$ok1 + $ok2;
 if ( $ok == 0 ) {
  open(TMP,">>OUT_SMFA");
  print TMP "calculation may have failed for $stem.mol\n";
  close TMP;
 }else{
  print DON "$stem.out\n";
 }
}
}











