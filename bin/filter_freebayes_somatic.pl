#!/usr/bin/perl -w
use strict;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use Data::Dumper;
use List::Util qw( max min );

my $vcf = vcf2->new('file'=>$ARGV[0] );

#TODO ADD INFO HEADER STRINGS PROPERLY!
system("zgrep ^## $ARGV[0]");
print '##INFO=<ID=SSC,Number=1,Type=Float,Description="Somatic score">'."\n";
system("zgrep ^#CHROM $ARGV[0]");

my $SSC_THRES = 30; # Somatic score threshold
my $LOD_THRES = 10; # Likelihood ratio threshold

while ( my $v = $vcf->next_var() ) {

    my @filters;


    my( %likelihood, %gl_idx, %genotype, %altobs, %depth );
    my $status = "PASS";
    for my $gt (@{$v->{GT}}) {
	my $type = (split /_/, $gt->{_sample_id})[-1];

	# Fail if GT is 0/0 for tumor
	$status = "FAIL_GT" if $type eq "T" and $gt->{GT} eq "0/0";
	
	my @GL = split /,/, $gt->{GL};
	my @GT = split /\//, $gt->{GT};
	my @AO = split /,/, ($gt->{AO} or "0");

	my $DP = $gt->{DP};
	my $RO = $gt->{RO};

	# This gets the position of the genotype's likelihood value in the GL array
	my $GL_IDX = ($GT[1]*($GT[1]+1)/2) + $GT[0];

	$depth{$type}      = $DP;
	$altobs{$type}     = \@AO;
	$genotype{$type}   = \@GT;
	$likelihood{$type} = \@GL;
	$gl_idx{$type} = $GL_IDX;
    }

    my $LOD_NORM  = $likelihood{N}->[$gl_idx{N}] - $likelihood{N}->[$gl_idx{T}];
    my $LOD_TUMOR = $likelihood{T}->[$gl_idx{T}] - $likelihood{T}->[$gl_idx{N}];
    my $DQUAL = $LOD_TUMOR+$LOD_NORM;

    # Loss of heterozygosity
    if( $genotype{T}->[0] eq $genotype{T}->[1] && ( $genotype{T}->[1] eq $genotype{N}->[0] or $genotype{T}->[1] eq $genotype{N}->[1] ) ) {
	$status = "LOH";
    }

    # Fail if low somatic score
    $status = "FAIL_QUAL" if $DQUAL<$SSC_THRES;

    # Fail if low likelyhood ratio hor either sample
    $status = "FAIL_LOD" if $LOD_NORM<$LOD_THRES || $LOD_TUMOR<$LOD_THRES;

    my $TALT;
    unless( $genotype{T}->[1] == $genotype{N}->[0] or $genotype{T}->[1] == $genotype{N}->[1] ) {
	$TALT = $genotype{T}->[1];
    }
    else {
	$TALT = $genotype{T}->[0];
	$status = "FAIL_NOVAR" if $TALT eq "0";
    }

    # Fail if difference between tumor's and normal's VAF is < 3x.
    if( $depth{N} > 0 and $depth{T} > 0 ) {
	my $NVAF = $altobs{N}->[$TALT-1] / $depth{N};
	my $TVAF = $altobs{T}->[$TALT-1] / $depth{T};
	if( $NVAF > 0 and ($TVAF/$NVAF < 3) ) {
	    $status = "FAIL_NVAF";
	}
    }

    $v->{FILTER} = $status;
    add_info( $v, "SSC", $DQUAL );    

    vcfstr($v);
}





sub add_info {
    my( $var, $key, $val ) = @_;
    push( @{$var->{INFO_order}}, $key );
    $var->{INFO}->{$key} = $val;
}

sub vcfstr {
    my $v = shift;

    my @all_info;
    print $v->{CHROM}."\t".$v->{POS}."\t".$v->{ID}."\t".$v->{REF}."\t".$v->{ALT}."\t".$v->{QUAL}."\t".$v->{FILTER}."\t";

    # Generate and print INFO field
    for my $info_key (@{$v->{INFO_order}}) {
	push @all_info, $info_key."=".$v->{INFO}->{$info_key};
    }
    print join(";", @all_info)."\t";

    # Print FORMAT field
    print join(":", @{$v->{FORMAT}})."\t";

    # Print GT fields for all samples
    for my $gt (@{$v->{GT}}) {
	my @all_gt;
	for my $key ( @{$v->{FORMAT}} ) {
	    push @all_gt, ( defined $gt->{$key} ? $gt->{$key} : "");
	}
	print join(":", @all_gt)."\t";
    }
    print "\n";
}
