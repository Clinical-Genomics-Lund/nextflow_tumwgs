#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;

my $MAX_PON_FRAC = 0.05;


# Get command line options
my %opt = ();
GetOptions( \%opt, 'vcf=s', 'tumor-id=s', 'pons=s' );
my %pon_files = check_options( \%opt );

# Get tumor ID
my $tid = $opt{'tumor-id'};

# Read all PONs into a hash
my %pons;
for my $vc (keys %pon_files) {
    $pons{$vc} = read_pon($pon_files{$vc});
}

# Open VCF
my $vcf = vcf2->new('file'=>$opt{vcf} );

# Print header
print_header($opt{vcf}, [keys %pons]);

while ( my $var = $vcf->next_var() ) {

    my $T_VAF;
    for my $gt ( @{$var->{GT}} ) {
	if( $tid and $gt->{_sample_id} eq $tid ) {
	    $T_VAF = $gt->{VAF};
	}
    }
    die "Tumor sample not found" if !defined($T_VAF);

    my $varid = $var->{CHROM}."_".$var->{POS}."_".$var->{REF}."_".$var->{ALT};

    my @filters = split ';', $var->{FILTER};
    my $remove_pass = 0;
    for my $vc (sort keys %pons) {
	my $pon = $pons{$vc};
	if( $pon->{$varid} ) {
	    my $p = $pon->{$varid};
	    my $frac = $p->{num_non_germ} / $p->{total};

	    if( $frac > $MAX_PON_FRAC ) {
		if( $T_VAF < $p->{mean_vaf} + $p->{stdev_vaf}*2 ) {
		    push @filters, "FAIL_PON_$vc";
		    $remove_pass = 1;
		}
		else {
		    push @filters, "WARN_PON_$vc";
		}
	    }
	    elsif( $p->{num_non_germ} >= 2 ) { # Add warning for anything that is in the PON at least twice
		push @filters, "WARN_PON_$vc";
	    }
	    add_info($var, "PON_NUM_$vc", $p->{num_non_germ}."/".$p->{total});
	    add_info($var, "PON_VAFS_$vc", $p->{all});
	
	}
    }

    @filters = grep {!/^PASS$/} @filters if $remove_pass;

    
    $var->{FILTER} = join(";", @filters);

    vcfstr($var);
}


sub read_pon {
    my $file = shift;
    
    open(PON, $file) or die "Cannot open $file";
    my %pon;
    while(<PON>) {
	chomp;
	my( $var, $vars, $germline, $total, $mean_vaf, $stdev_vaf, $all_vafs) = split /\t/;
	$pon{$var} = {
	    'num_non_germ'=>($vars-$germline),
	    'total'=>$total,
	    'mean_vaf'=>$mean_vaf,
	    'stdev_vaf'=>$stdev_vaf,
	    'all'=>$all_vafs
	};
}

    return \%pon;
}

sub print_header {
    my $file = shift;
    my $vcs = shift;

    system("zgrep ^## $file");
    for my $vc (@$vcs) {
	print "##FILTER=<ID=WARN_PON_$vc,Description=\"Warning, exists in PON for $vc\">\n";
	print "##FILTER=<ID=FAIL_PON_$vc,Description=\"Failed, exists in PON for $vc\">\n";

	print "##INFO=<ID=PON_NUM_$vc,Number=.,Type=String,Description=\"Number of PON samples variants exists in, for $vc\">\n";
	print "##INFO=<ID=PON_VAFS_$vc,Number=.,Type=String,Description=\"VAFs for PON variants, for $vc\">\n";
    }
   
    system("zgrep ^#CHROM $file");
	    
}




sub check_options {
    my %op = %{ $_[0] };

    help_text() unless $op{vcf} and $op{pons};

    my %pons;
    my @pons = split /,/, $op{pons};
    foreach( @pons ) {
	my( $vc, $file) = split /=/;
	$pons{$vc} = $file;
	die "PON file does not exist $file..." unless -s $file;
    }
    
    die "File does not exist $op{vcf}..." if ! -s $op{vcf};

    return %pons;

}


sub help_text {
    my $error = shift;
    
    print "\n\$ mark_germlines.pl --vcf INPUT_VCF --tumor-id [--normal-id}\n\n";
    print "   --vcf        Input vcf\n";
    print "   --tumor-id   Tumor sample ID\n";
    print "   --pon        Path to PON file, created with create_snv_pon.pl\n";
    print "\n";
    exit(0);
}


sub vcfstr {
    my( $v ) = @_;

    my @all_info;
    print $v->{CHROM}."\t".$v->{POS}."\t".$v->{ID}."\t".$v->{REF}."\t".$v->{ALT}."\t".$v->{QUAL}."\t".$v->{FILTER}."\t";

    # Generate and print INFO field
    for my $info_key (@{$v->{INFO_order}}) {
	my $key2 = $info_key;
	$key2 = "_CSQ_str" if $info_key eq "CSQ";
	push @all_info, $info_key."=".$v->{INFO}->{$key2};
    }
    print join(";", @all_info)."\t";

    # Print FORMAT field
    print join(":", @{$v->{FORMAT}});

    
    # Print GT fields for all samples
    for my $gt (@{$v->{GT}}) {
	my @all_gt;
	for my $key ( @{$v->{FORMAT}} ) {
	    push @all_gt, ( defined $gt->{$key} ? $gt->{$key} : "");
	}
	print "\t".join(":", @all_gt);
    }
    print "\n";
}

sub add_info {
    my( $var, $key, $val ) = @_;
    push( @{$var->{INFO_order}}, $key );
    $var->{INFO}->{$key} = $val;
}
