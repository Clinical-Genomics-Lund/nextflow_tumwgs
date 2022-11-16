#!/usr/bin/perl -w
use strict;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use Data::Dumper;

#/fs1/results_dev/tumwgs/vcf/7992-19.agg.pon.vep.vcf
my $vcf_fn = $ARGV[0];
my $panel_fn = $ARGV[1];
my $hard_filter = $ARGV[2] ? 1 : 0;

my %panel = read_panel($panel_fn);

my $vcf = vcf2->new('file'=>$ARGV[0] );

#TODO ADD INFO HEADER STRINGS PROPERLY!
system("zgrep ^## $vcf_fn");
system("zgrep ^#CHROM $vcf_fn");

while ( my $v = $vcf->next_var() ) {
    if (is_in_panel($v, \%panel) or !$hard_filter) {
	    vcfstr($v);
    }
    else {
         vcfstr($v);   
    }
}


sub is_in_panel {
    my( $v, $panel ) = @_;

    my %classes;
    foreach my $csq ( @{$v->{INFO}->{CSQ}} ) {
	my $gene = $csq->{SYMBOL};
	next if $csq->{IMPACT} eq "MODIFIER";
	if( $panel->{$gene} ) {
	    $classes{$gene} = $panel->{$gene};
	}
    }

    if( %classes ) {
	add_info($v, 'PANEL', join(",",values %classes));
	return 1;
    }

    return 0;
}

sub read_panel {
    my $fn = shift;

    my %panel;
    open( PANEL, $fn );
    while(<PANEL>) {
	chomp;
	my( $gene, $class, $snv ) = split /\t/;
	$panel{$gene} = $class.":".($snv == 1 ? "snv" : "other");
    }
    return %panel;
}
       
sub add_info {
    my( $var, $key, $val ) = @_;
    push( @{$var->{INFO_order}}, $key );
    $var->{INFO}->{$key} = $val;
}


sub vcfstr {
    my( $v, $sample_order ) = @_;

    my @all_info;
    print $v->{CHROM}."\t".$v->{POS}."\t".$v->{ID}."\t".$v->{REF}."\t".$v->{ALT}."\t".$v->{QUAL}."\t".$v->{FILTER}."\t";

    # Generate and print INFO field
    for my $info_key (@{$v->{INFO_order}}) {
	    if( $info_key eq "CSQ" ) {
	        push @all_info, "CSQ=".$v->{_csqstr};
	    }
	    else {
	        push @all_info, $info_key."=".$v->{INFO}->{$info_key};
	    }
    }
    print join(";", @all_info)."\t";

    # Print FORMAT field
    print join(":", @{$v->{FORMAT}})."\t";


    my %order;
    my $i=0;
    if( $sample_order and @$sample_order > 0 ) {
	    $order{$_} = $i++ foreach @{$sample_order};
    }
    else {
	    $order{$_->{_sample_id}} = $i++ foreach @{$v->{GT}};
    }

    # Print GT fields for all samples
    for my $gt ( sort {$order{$a->{_sample_id}} <=> $order{$b->{_sample_id}}} @{$v->{GT}}) {
	
    my @all_gt;
	for my $key ( @{$v->{FORMAT}} ) {
	    push @all_gt, ( defined $gt->{$key} ? $gt->{$key} : "");
	}
	print join(":", @all_gt)."\t";
    }
    print "\n";
}
