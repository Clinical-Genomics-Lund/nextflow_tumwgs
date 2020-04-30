#!/usr/bin/perl -w
use strict;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use Data::Dumper;

#/fs1/results_dev/tumwgs/vcf/7992-19.agg.pon.vep.vcf
my $vcf_fn = $ARGV[0];
my $panel_fn = $ARGV[1];

my %panel = read_panel($panel_fn);

my $vcf = vcf2->new('file'=>$ARGV[0] );

#TODO ADD INFO HEADER STRINGS PROPERLY!
system("zgrep ^## $vcf_fn");
system("zgrep ^#CHROM $vcf_fn");

my %partner_seen;
while ( my $v = $vcf->next_var() ) {
    next unless $v->{INFO}->{SVTYPE} eq "BND";

    next if $partner_seen{$v->{CHROM}.":".$v->{POS}};
    
    my( $partner_chr, $partner_pos ) = parse_partner($v->{ALT});
    $partner_seen{$partner_chr.":".$partner_pos} = 1;
    
    my $max_class = 0;
    my $in_both = 0;
    foreach my $ann ( @{ $v->{INFO}->{ANN} } ) {
	my @genes = split /&/, $ann->{Gene_Name};
	next if @genes < 2;

	my( $class, $both ) = is_in_panel($v, \%panel, $genes[0], $genes[1]);

	$in_both = 1 if $both;
	$max_class = max_class($class, $max_class);
    }

    if( $max_class ) {
	add_info($v, 'PANEL', 'fusion|'.$max_class.'|'.($in_both ? 'both' : 'one'));
	vcfstr($v);
    }
}

sub max_class {
    my( $c1, $c2 ) = @_;
    
    return "somatic" if $c1 eq "somatic" or $c2 eq "somatic";
    return "research" if $c1 eq "research" or $c2 eq "research";
    return "none" if $c1 eq "none" or $c2 eq "none";
    return 0;
}

sub is_in_panel {
    my( $v, $panel, $g1, $g2 ) = @_;

    my $class = 0;
    my $both = 0;
    
    if( $panel->{$g1} ) {
	if( defined $panel->{$g1}->{$g2} ) {
	    $class = ( $panel->{$g1}->{$g2} or "none" );
	    $both = 1;
	}
	else {
	    $class = "none";
	}
    }
    elsif( $panel->{$g2} ) {
	if( defined $panel->{$g2}->{$g1} ) {
	    $class = ( $panel->{$g2}->{$g1} or "none" );
	    $both = 1;
	}
	else {
	    $class = "none";
	}
    }

    return($class, $both);
    
}

sub read_panel {
    my $fn = shift;

    my %panel;
    open( PANEL, $fn );
    while(<PANEL>) {
	chomp;
	my( $gene, $partners, $class ) = split /\t/;
	my @partners = split /,/, $partners;
	foreach my $partner ( @partners ) {
	    $panel{$gene}->{$partner} = $class;
	}
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
	if( $info_key eq "ANN" ) {
	    push @all_info, "ANN=".$v->{_annstr};
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

sub parse_partner {
    my $alt = shift;

    if( $alt =~ /\[(.*?):(.*?)\[/ ) {
	return( $1, $2 );
    }
    if( $alt =~ /\](.*?):(.*?)\]/ ) {
	return( $1, $2);
    }
    return(0,0);
}
