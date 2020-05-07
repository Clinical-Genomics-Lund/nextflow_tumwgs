#!/usr/bin/perl -w
use strict;

my $in_bed = $ARGV[0];
my $genes_bed = $ARGV[1];

my $rnd = int rand 1000000000;
my $tmp_infile = "input.$rnd.bed";
system("grep -v '^\@' $in_bed| grep -v ^CONTIG > $tmp_infile");
my @overlap = `bedtools intersect -a $tmp_infile -b $genes_bed -loj`;
unlink $tmp_infile;

my $prev_reg;
my @genes;
my $prev_bed_str;
foreach my $line (@overlap) {
    chomp $line;
    
    my @f = split /\t/, $line;

    my $reg = "$f[0]\t$f[1]\t$f[2]\t";
    my $bed_str = "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]";
    
    if( $prev_reg and $reg ne $prev_reg ) {
        print $prev_bed_str."\t".join(",", @genes)."\n";
        @genes = ();
    }
    push @genes, $f[-1];
    
    $prev_reg = $reg;
    $prev_bed_str = $bed_str;
}

print $prev_bed_str."\t".join(",", @genes)."\n";

#    chr116154141692795chr116154141630610MIB2
#    chr116154141692795chr116352261659012CDK11B
#    chr116154141692795chr116595281692795SLC35E2B
#    chr210001002.-1-1.
#    chr1016154142692795chr1019575881960190AL355597.1
#    chr1016154142692795chr1020760182076089MIR6072
#    chr1016154142692795chr1020981662098269RNU6-889P
#    chr1016154142692795chr1021993862199488RNU6-576P
