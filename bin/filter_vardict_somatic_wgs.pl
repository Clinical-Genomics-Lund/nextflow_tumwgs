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

my $MIN_VAF_RATIO = 3;
my $MIN_VAF_HOMOPOLYMER_RATIO = 5;

#TODO ADD INFO HEADER STRINGS PROPERLY!
system("zgrep ^## $ARGV[0]");
system("zgrep ^#CHROM $ARGV[0]");


while ( my $v = $vcf->next_var() ) {

    #    my $status = "PASS";
    my @status;
    my (%vaf, %dp, %strand_bias, %msi, %msilen);

    my $msi    = $v->{INFO}->{MSI};
    my $msilen = $v->{INFO}->{MSILEN};

    my $is_indel = 0;
    $is_indel = 1 if length($v->{REF}) != length($v->{ALT});
    
    for my $gt (@{$v->{GT}}) {
        my $type = "T";
        $type = "N" if $gt->{_sample_id} eq $N;

	$vaf{$type} = $gt->{AF};
	$dp{$type} = $gt->{DP};
	$strand_bias{$type} = $gt->{SBF};
    }

    # Fail long INDELs as they tend to be false positives in VarDict.
    if( length($v->{REF}) - length($v->{ALT}) > 250 or $v->{ALT} =~ /^<DEL>$/ ) {
	push @status, "FAIL_LONGDEL";
    }
    elsif( length($v->{ALT}) - length($v->{REF}) > 250 or $v->{ALT} =~ /^<(DUP|INS)>$/ ) {
	push @status, "FAIL_LONGINS";
    }
    
    # Fail if difference between tumor's and normal's VAF is < 3x.
    if( $vaf{T} > 0 ) {
	if( $vaf{N} > 0 and ($vaf{T} / $vaf{N} < $MIN_VAF_RATIO) ) {
	    push @status, "FAIL_NVAF";
	}
	if( $dp{T} < 30 ) {
	    push @status, "WARN_LOW_TCOV";
	}
	if( $vaf{T} < 0.02 ) {
	    push @status, "FAIL_VERYLOW_TVAF";
	}
	elsif( $vaf{T} < 0.05 ) {
	    push @status, "FAIL_LOW_TVAF";
	}
	if( $strand_bias{T} < 0.05 ) {
	    push @status, "WARN_STRANDBIAS";
	}

        # If in a homopolymer
	if( $msilen <= 2 and $msi > 10 ) {
	    if( !$is_indel or !$vaf{N} or $vaf{T} / $vaf{N} >= $MIN_VAF_HOMOPOLYMER_RATIO ) {
		push @status, "WARN_HOMOPOLYMER_".($is_indel?"INDEL":"SNV");
	    }
	    else {
		push @status, "FAIL_HOMOPOLYMER_".($is_indel?"INDEL":"SNV");
	    }
	}

	
	
    }
    else {
	push @status, "FAIL_NO_TVAR";
    }

    if( @status ) {
	if( $v->{FILTER} eq "PASS" or $v->{FILTER} eq "." ) {
	    $v->{FILTER} = join(";", @status);
	}
	else {
	    $v->{FILTER} .= ";".join(";", @status);
	}
    }
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
