#!/usr/bin/env perl
use v5.14;
use warnings;
use strict;
my $split_len=shift;
my $size = $split_len;
my $j=0;
my $i=1;
while(<>){
    chomp;
    my $id = $_ if(/>/);
    my $seq = <>;
    $j+=length($seq);
    if($j<=$size){
	open(O,">>Split-$i.fa");
	print O "$id\n$seq";
	close O;
    }else{
	$i++;$j=length($seq);
	open(O,">>Split-$i.fa");
	print O "$id\n$seq";close O;
    }
}
