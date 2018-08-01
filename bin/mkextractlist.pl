#!/usr/bin/perl

$CODEDIR=`awk 'NR==1,NR==1 {print \$1}' CODENAME`;
chomp($CODEDIR);
$CODEDIR=~ s/^\s+|\s+$//g;

$filename='ABSETTINGS';
open($fh,'<:encoding(UTF-8)',$filename)
or die "could not open '$filename' $!";
$n=0;
while ($row = <$fh>){
chomp $row;
$store[$n]=$row;
$n++;
}
$n=$n-1;
close $fh;

$qmprog=$store[1];

$nfrgs = `ls COORD* | wc -l`;
chomp($nfrgs);
$nfrgs=~ s/^\s+|\s+$//g;

$line=`awk /ncpus=/ pbsfile`;
@bline=split("=",$line);
$ncpus=$bline[1];
chomp($ncpus);
$ncpus=~ s/^\s+|\s+$//g;

$div=int($nfrgs / $ncpus);
$rem= $nfrgs -  $div * $ncpus;
@llim[1]=1;
if($rem == 0) {
 @ulim[1]=$div;
}else{
 @ulim[1]=$div+1;
}
for ($i=2;$i<=$ncpus;$i++) {
 $j=$i-1;
 if($i <= $rem) {
  @llim[$i]=@ulim[$j]+1;
  @ulim[$i]=@llim[$i]+$div;
 }else{
  @llim[$i]=@ulim[$j]+1;
  @ulim[$i]=@llim[$i]+$div-1;
 }
}

if($qmprog == 1) {$scr="extract_gam.pl"};
if($qmprog == 2) {$scr="extract_gau.pl"};
if($qmprog == 3) {$scr="extract_nwc.pl"};
if($qmprog == 4) {$scr="extract_qch.pl"};

$junk=`rm -f getdata*`;

open(OUT,">getdata");
for ($i=1;$i<=$ncpus;$i++) {
print OUT "$CODEDIR/$scr  $llim[$i]    $ulim[$i]   $i  &\n";
}
print OUT "wait\n";
close OUT;

&multnodeoutput("getdata");
system("chmod +x getdata");

# do the ab file lists

$junk=`ls -1 ab.*.log > abfiles`;

$line=`awk /ncpus=/ pbsfile`;
@bline=split("=",$line);
$ncpus=$bline[1];
chomp($ncpus);
$ncpus=~ s/^\s+|\s+$//g;

$junk=`rm -f abfiles.*`;
$ic=0;
open(IN,"<abfiles");
while (<IN>) {
$ic=$ic+1;
$line=$_;
chomp($line);
$line=~ s/^\s+|\s+$//g;
if ($ic > $ncpus) {$ic = $ic - $ncpus};
system("echo $line >> abfiles.$ic");
}


sub multnodeoutput {
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

$junk=`rmdir nbdr* 2> /dev/null`;


