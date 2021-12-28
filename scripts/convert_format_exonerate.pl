#/usr/bin/perl -w
#perl ~/bin/annation/exonerate2evm.pl  Split-2.fa_out  | awk  'BEGIN{i=0};{if($0~/alignment_id 1\s/){i++};print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t" "ID=match_NO_" i ";AveragePercentIdentity="$19 }' | grep -v "line" 
open IN ,shift;
my $num=0;
foreach my $file (<IN>){
    chomp $file ;
    if ($file =~/alignment_id/){
        my @file = split (/\t/,$file);
        print "$file[0]\t$file[1]\tnucleotide_to_protein_match\t$file[3]\t$file[4]\t$file[5]\t$file[6]\t$file[7]\t$file[8];"
    }
    if($file =~/AveragePercentIdentity:/){
        print "$file\n";
    }
}
