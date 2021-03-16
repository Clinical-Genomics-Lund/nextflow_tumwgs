#!/usr/bin/perl -w
use strict;
use File::Basename;
use lib dirname (__FILE__);
use Data::Dumper;

my $bed_fn = $ARGV[0];
my $panel_fn = $ARGV[1];

my %panel = read_panel($panel_fn);

open(BED, $bed_fn);
while ( my $v = <BED> ) {
    chomp $v;
    my @v = split /\t/, $v;
    next if $v[4] > -0.25 and $v[4] < 0.25;

    my $panel_info = in_panel($v[6], \%panel);
    if( $panel_info) {
	print $v."\t".$panel_info."\n";
    }
    
}
close BED;

sub in_panel {
    my( $genes, $panel ) = @_;

    my @genes = split /,/, $genes;
    my %classes;
    foreach my $gene ( @genes ) {
	if( $panel->{$gene} ) {
	    $classes{$gene} = $gene.":".$panel->{$gene};
	}
    }

    if( %classes ) {
	return join(",",values %classes);
    }

    return 0;
}

sub read_panel {
    my $fn = shift;

    my %panel;
    open( PANEL, $fn );
    while(<PANEL>) {
	chomp;
	my( $gene, $type, $class ) = split /\t/;
	$panel{$gene} = $class.":".$type;
    }
    return %panel;
}
