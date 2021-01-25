#!/usr/bin/perl -w
use strict;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use Data::Dumper;
use List::Util qw( max min );

my $vcf = vcf2->new('file'=>$ARGV[0] );

my $T = $ARGV[1];
my $N = $ARGV[2];

#TODO ADD INFO HEADER STRINGS PROPERLY!
system("zgrep ^## $ARGV[0]");
print '##INFO=<ID=SSC,Number=1,Type=Float,Description="Somatic score">'."\n";
system("zgrep ^#CHROM $ARGV[0]");

my $SSC_THRES = 10; # Somatic score threshold
my $LOD_THRES = 5; # Likelihood ratio threshold
my $MIN_VAF_RATIO = 3;

while ( my $v = $vcf->next_var() ) {

    my @filters;


    my( %likelihood, %gl_idx, %genotype, %altobs, %depth, %refobs, %alt_qualsum, %ref_qualsum);
    my $status = "PASS";
    for my $gt (@{$v->{GT}}) {
        my $type = "T";
        $type = "N" if $gt->{_sample_id} eq $N;

	# Fail if GT is 0/0 for tumor
	if( $gt->{GT} =~ /\./ ) {
	    $status = "FAIL_GT";
	    next;
	}
	
	my @GL = split /,/, $gt->{GL};
	my @GT = split /\//, $gt->{GT};
	my @AO = split /,/, ($gt->{AO} or "0");

	my $DP = $gt->{DP};
	my $RO = $gt->{RO};
	my @QA = split /,/, ($gt->{QA} or "0");
	my $QR = $gt->{QR};

	# This gets the position of the genotype's likelihood value in the GL array
	my $GL_IDX = ($GT[1]*($GT[1]+1)/2) + $GT[0];

	$depth{$type}      = $DP;
	$altobs{$type}     = \@AO;
	$refobs{$type}     = $RO;
	$genotype{$type}   = \@GT;
	$likelihood{$type} = \@GL;
	$alt_qualsum{$type} = \@QA;
	$ref_qualsum{$type} = $QR;
	$gl_idx{$type} = $GL_IDX;
    }

    unless( $status eq "FAIL_GT" ) {
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
	$status = "FAIL_LOD.$LOD_TUMOR.$LOD_NORM" if $LOD_NORM<$LOD_THRES || $LOD_TUMOR<$LOD_THRES;
	
	my $TALT;
	unless( $genotype{T}->[1] == $genotype{N}->[0] or $genotype{T}->[1] == $genotype{N}->[1] ) {
	    $TALT = $genotype{T}->[1];
	}
	else {
	    $TALT = $genotype{T}->[0];
	    $status = "WARN_NOVAR" if $TALT eq "0";
	}
	

	# Fail if average base qual score is much lower for ALT than for REF 
	if( $ref_qualsum{T} ) {
	    my $avg_qual_ref = $ref_qualsum{T}/$refobs{T};
	    my $avg_qual_alt = $alt_qualsum{T}->[$TALT-1] / $altobs{T}->[$TALT-1];
	    my $refalt_qratio = $avg_qual_ref/$avg_qual_alt;
	    add_info( $v, "RA_Qratio", $refalt_qratio);
	    $status = "FAIL_LOWALTQ" if $refalt_qratio >= 2.5;
	    $status = "FAIL_LOWALTQ" if $refalt_qratio >= 1.8 and $avg_qual_alt < 17;
	}
	
	
	# Fail if difference between tumor's and normal's VAF is < 3x.
	if( $depth{N} > 0 and $depth{T} > 0 ) {
	    my $NVAF = $altobs{N}->[$TALT-1] / $depth{N};
	    my $TVAF = $altobs{T}->[$TALT-1] / $depth{T};

	    if( $TVAF < 0.05 ) {
		$status = "FAIL_LOW_TVAF";
	    }
	    if( $NVAF > 0 and ($TVAF/$NVAF < $MIN_VAF_RATIO) ) {
		$status = "FAIL_NVAF";
	    }
	    
	}

	# Warn if variant was only found on one strand
	if( $status !~ /^FAIL/) {
	    my @SAF = split /,/, $v->{INFO}->{SAF};
	    my @SAR = split /,/, $v->{INFO}->{SAR};
	    $status = "WARN_STRAND.$status" if $SAF[$TALT-1] == 0 or $SAR[$TALT-1] == 0;
	}
	

	
	add_info( $v, "SSC", $DQUAL );


    }

    $v->{FILTER} = $status;

    
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
