#!/usr/bin/perl

$list=$ARGV[0];
chomp($list);
$list=~ s/^\s+|\s+$//g;
$listorig=join("",$list,"_original");
system("cp $listorig $list");

@a=split("IN",$list);
$a1=$a[1];
$out=join("","OUT",$a1);

open(OUT,"<$out");
while (<OUT>) {
$line=$_;
@a=split("log",$line);
$stem=$a[0];
system("sed -i '/$stem/d' ./$list");
}


