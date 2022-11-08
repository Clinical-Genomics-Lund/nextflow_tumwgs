#!/usr/bin/env perl
use strict;
use warnings;
use vcf2;
use Data::Dumper;

my $TID = "10210-18";
my $NID = "10213-18";
my $in = "/fs1/results_dev/tumwgs/manta/277-20_manta.vcf.snpeff";

my $vcf = vcf2->new('file'=>$in );

my $cnt;
my %bnds;
while ( my $v = $vcf->next_var() ) {
    if( $v->{INFO}->{SVTYPE} eq "BND") {
	my %pr_sr = parse_PR_SR($v->{GT});
	my( $partner_chr, $partner_pos ) = parse_partner($v->{ALT});
	
	my $bnd_id = ($v->{INFO}->{EVENT} or join( "_", sort($v->{CHROM}.":".$v->{POS}, $partner_chr.":".$partner_pos)));
	my $pos = $v->{CHROM}.":".$v->{POS};
#	unless( $bnds{$bnd_id} ) {
#	    $bnds{$bnd_id}->{POS} = $pos;
#	}
	push @{$bnds{$bnd_id}->{POS}}, $pos;

	my @consequences = parse_consequences($v->{INFO}->{ANN});
	$bnds{$bnd_id}->{conseq} = \@consequences;

	$bnds{$bnd_id}->{hgvsp} = collect_from_ann($v->{INFO}->{ANN}, "HGVS.p");
	$bnds{$bnd_id}->{genes} = collect_from_ann($v->{INFO}->{ANN}, "Gene_Name");

	#$v->{INFO}->{CSQ} = 0;
	#print Dumper($v);
	$cnt++;
    }
}

print Dumper(\%bnds);

sub parse_consequences {
    my $anns = shift;

    my %csqs;
    foreach my $ann (@$anns) {
	foreach my $consequence (@{$ann->{Annotation}}) {
	    $csqs{$consequence}++;
	}
    }
    return (keys %csqs);
}

sub collect_from_ann {
    my $anns = shift;
    my $field = shift;

    my %vals;
    foreach my $ann (@$anns) {
	$vals{$ann->{$field}}++;
    }
    return [keys %vals];
}

sub parse_PR_SR {
    my $gts = shift;

    my %data;
    foreach my $gt (@$gts) {
	for my $type ("SR", "PR") {
	    if( $gt->{$type} ) {
		my( $ref,$alt ) = split /,/, $gt->{$type};
		$data{$gt->{_sample_id}}->{$type} = {'ref'=>$ref, 'alt'=>$alt};
	    }
	}
    }
    return %data;
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
    
