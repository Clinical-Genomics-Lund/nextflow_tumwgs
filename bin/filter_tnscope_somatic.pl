#!/usr/bin/perl -w
use strict;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use Data::Dumper;
use List::Util qw( max min sum );

#1       1718690 .       A       T       0.14    ML_FAIL ECNT=1;FS=0;HCNT=8;MAX_ED=.;MIN_ED=.;ML_PROB=0.02;NLOD=181.64;NLODF=54.2;PV=0.085;PV2=0.0561;SOR=0.415;TLOD=4.81        GT:AD:AF:AFDP:ALTHC:ALT_F1R2:ALT_F2R1:BaseQRankSumPS:ClippingRankSumPS:DPHC:FOXOG:MQRankSumPS:NBQPS:QSS:REF_F1R2:REF_F2R1:ReadPosEndDistPS:ReadPosRankSumPS     0/1:715,14:0.02:703:12:12:2:-7.451:0:720:0.143:0:25.963:21370,330:368:347:30.542:-5.755 0/0:696,6:0.009:646:4:5:1:-4.974:0:629:0.167:0:26.017:20615,119:354:342:30.751:-2.414
   

my $vcf = vcf2->new('file'=>$ARGV[0] );

my $T = $ARGV[1];
my $N = $ARGV[2];

my $MIN_VAF_RATIO = 3;
my $MIN_VAF_HOMOPOLYMER_RATIO = 5;

#TODO ADD INFO HEADER STRINGS PROPERLY!
system("zgrep ^## $ARGV[0]");
system("zgrep ^#CHROM $ARGV[0]");


while ( my $v = $vcf->next_var() ) {

    my @status;
    my (%vaf, %dp);

    my $is_indel = 0;
    $is_indel = 1 if length($v->{REF}) != length($v->{ALT});

    my( $repeat_units_num, $repeat_seq ) = (0, "");
    if( $v->{INFO}->{RPA} and $v->{INFO}->{RU} ) {
	$repeat_units_num = min( split /,/, $v->{INFO}->{RPA} );
	$repeat_seq = ($v->{INFO}->{RU} or "");
    }
    	    
    
    for my $gt (@{$v->{GT}}) {
        my $type = "T";
        $type = "N" if $gt->{_sample_id} eq $N;

	$vaf{$type} = $gt->{AF};
	$dp{$type} = sum( split /,/, $gt->{AD} );
	    
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
	    push @status, "WARN_VERYLOW_TVAF";
	}
	elsif( $vaf{T} < 0.05 ) {
	    push @status, "WARN_LOW_TVAF";
	}

	# If in a homopolymer or dinucleotide repeat
	if( length($repeat_seq) <= 2 and $repeat_units_num >= 10 ) {
	    if( $is_indel or !$vaf{N} or $vaf{T} / $vaf{N} >= $MIN_VAF_HOMOPOLYMER_RATIO ) {
		push @status, "WARN_HOMOPOLYMER_".($is_indel?"INDEL":"SNV");
	    }
	    else {
		push @status, "FAIL_HOMOPOLYMER_".($is_indel?"INDEL":"SNV");
	    }
	}

	# Warn/fail on high strand bias (SOR)
	if( $is_indel and $v->{INFO}->{SOR} > 4) {
	    push @status, ( $v->{INFO}->{SOR} > 10 ? "FAIL_STRANDBIAS" : "WARN_STRANDBIAS" );
	}
	if( !$is_indel and $v->{INFO}->{SOR} > 2.5) {
	    push @status, ( $v->{INFO}->{SOR} > 4 ? "FAIL_STRANDBIAS" : "WARN_STRANDBIAS" );
	}

	if( $v->{INFO}->{PV2} > 0.05 ) {
	    push @status, "FAIL_PVALUE";
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
