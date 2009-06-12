#!/usr/bin/perl -w

# Converts a CAF file into a BAF file. Accepts padded data only

use strict;
use Caftools;

#---------- Load CAF file
my $caf = Caftools->new();
$caf->readCafFile(\*STDIN);

#---------- Iterate through readings in contigs producing BAF
foreach my $contig ($caf->contigsList) {
    local $"="";

    my $entry = $caf->get_entry_ref($contig);
    my %s;
    foreach (split("\n",${$entry->{Sequence}})) {
	my ($k,$v) = /(\S+)\s*(.*)/;
	$v =~ s/^"(.*)"$/$1/;
	if ($k ne "Align_to_SCF" && $k ne "Tag") {
	    $s{$k}=$v;
	} else {
	    push(@{$s{$k}}, $v);
	}
    }

    print "CO=$contig\n";
    print "LN=", length(${$entry->{DNA}}), "\n";
    print "\n";

    # Consensus tags
    if (exists($s{Tag})) {
	foreach (@{$s{Tag}}) {
	    my ($type,$st,$en,$text) =
		m/(\S+)\s+(\d+)\s+(\d+)\s*(?:"(.*)")?/;
	    print "AN=$type\n";
	    if ($en >= $st) {
		print "LO=$st\n";
		print "LL=",$en-$st+1,"\n";
	    } else {
		print "LO=$en\n";
		print "LL=",$st-$en+1,"\n";
	    }
	    print "TX=$text\n" if defined($text);
	    print "\n";
	}
    }
    

    #------ Reads within the contigs
    foreach my $af (@{$caf->get_assembled_from($contig)}) {

	# Data from the AssembledFrom records
	my ($read,$cstart,$cend,$rstart,$rend) = @$af;
	my $dir = $cstart < $cend ? 1 : -1;
	if ($dir == -1) {
	    my $tmp = $cstart;
	    $cstart = $cend;
	    $cend = $tmp;
	}
	print "RD=$read\n";
	print "DR=$dir\n";
	print "AP=$cstart\n";
	print "QL=$rstart\n";
	print "QR=$rend\n";

	# Data from "Sequence:" blocks
	$entry = $caf->get_entry_ref($read);
	my %s;
	foreach (split("\n",${$entry->{Sequence}})) {
	    my ($k,$v) = /(\S+)\s*(.*)/;
	    $v =~ s/^"(.*)"$/$1/;
	    if ($k ne "Align_to_SCF" && $k ne "Tag") {
		$s{$k}=$v;
	    } else {
		push(@{$s{$k}}, $v);
	    }
	}
	if (exists($s{Strand})) {
	    print "PR=0\n" if $s{Strand} eq "Forward";
	    print "PR=1\n" if $s{Strand} eq "Reverse";
	}
	print "TN=$s{Template}\n" if exists($s{Template});
	print "TR=$s{SCF_File}\n" if exists($s{SCF_File});

	die "Please run this file through caf_pad first"
	    unless exists $s{Padded};

	# Alignment edit string from the Align_to_SCF records.
	# This is a bit like a cigar string with M, I and D elements, but
	# we also have an E element marking "edit". In a compound case
	# we cannot distinguish by edit followed by indel and indel followed
	# by edit, so we guess.
	my $al = "";
	my $gaps = 0;
	my $last_sr = 0;
	my $last_tr = 0;
	if (1) {
	    foreach (@{$s{Align_to_SCF}}) {
		my ($sl,$sr,$tl,$tr) = split(" ", $_);

		my $ds = $sl-$last_sr-1;
		my $dt = $tl-$last_tr-1;

		if ($ds > $dt) {
		    if ($dt) {
			$al .= al_str("E", $dt);
			$ds -= $dt;
		    }
		    $al .= al_str("I", $ds);
		} elsif ($ds < $dt) {
		    if ($ds) {
			$al .= al_str("E", $ds);
			$dt -= $ds;
		    }
		    $al .= al_str("D", $dt);
		} else {
		    if ($ds) {
			$al .= al_str("E", $ds);
		    }
		}
		$al .= al_str("M", $sr-$sl+1);

		$last_sr = $sr;
		$last_tr = $tr;
	    }
	} else {
	    # A more compact format:
	    #
	    # Align_to_SCF can be compressed to size of range, gap in
	    # sequence coordinates, gap in trace coorindates. This triple is
	    # sufficient to encode all records.
	    my $last_sr = 0;
	    my $last_tr = 0;
	    foreach (@{$s{Align_to_SCF}}) {
		my ($sl,$sr,$tl,$tr) = split(" ", $_);

		my $sz = $sr-$sl+1;
		my $ds = $sl-$last_sr-1;
		my $dt = $tl-$last_tr-1;
		#$al .= "$sz,$ds,$dt";
		$al .= al_val($sz,$ds,$dt);
		
		$last_sr = $sr;
		$last_tr = $tr;
	    }
	}
	print "AL=$al\n";

	# Seq and quality
	print "SQ=${$entry->{DNA}}\n";
	my @bq = map {chr(33+$_)} split(" ", ${$entry->{BaseQuality}});
	print "FQ=@bq\n";

	print "\n";

	# Annotations
	if (exists($s{Tag})) {
	    foreach (@{$s{Tag}}) {
		my ($type,$st,$en,$text) =
		    m/(\S+)\s+(\d+)\s+(\d+)\s*(?:"(.*)")?/;
		print "AN=$type\n";
		if ($en >= $st) {
		    print "LO=$st\n";
		    print "LL=",$en-$st+1,"\n";
		} else {
		    print "LO=$en\n";
		    print "LL=",$st-$en+1,"\n";
		}
		print "TX=$text\n" if defined($text);
		print "\n";
	    }
	}
    }
}

sub al_str {
    my ($type, $size) = @_;

    return $type . $size;

    # Using a larger base to shrink numbers saves 10-15% of AL size,
    # but is largely irrelevant overall.
    my $b64 = "";
    while ($size > 0) {
	$b64 .= chr(33+ ($size%64));
	$size = int($size/64);
    }
    return lc($type) . $b64;
}

# ASCII-fy a number.
# 5 bits for values, top bit for continuation (more to come). All with
# +33 to make printable
sub al_val {
    my $v = "";
    foreach my $sz (@_) {
	do {
	    if ($sz >= 32) {
		$v .= chr(33+($sz%32)+32);
	    } else {
		$v .= chr(33+$sz);
	    }
	    $sz = int($sz/32);
	} while ($sz > 0);
    }
    return $v;
}
