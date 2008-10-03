#!/usr/bin/perl -w

# This attempts to reorder an ACE file so that the reads within it are
# sorted in left to right order.
#
# The purpose for this is to improve efficiency when converting to a Gap5
# database. While it'll work regardless of the read order (provided the
# reads are presented in the same order as listed int he AF lines), it is
# inefficient if the reads for a small region are spread throughout a large
# area of disk.

use strict;

# work in progress. Rewrite phrap2caf as it's SO slow as to be almost
# unusable.

my %seq;
my @AF;
my @BS;

my $curr_seq = "";

my $cnum = 1;
while (<>) {
    if (/^AS /) {
	print;
	$_=<>; # swallow blank line too
	print;
    }

    # RD D8QSF4J01CU6TZ 96 0 0
    # AGGGAGGG*AAG*AAATACACCTATG*CTTCCTCAGTAGGTAAAG*****
    # ***AAGA*C*TGGG***AACC*ACCACG*CCAGAGTTA*GAAAATA
    if (/^RD (\S+)/) {
	$curr_seq = $1;
	$seq{$curr_seq} = $_;
	while (<>) {
	    $seq{$curr_seq} .= $_;
	    last if (/^$/);
	}
    }

    # CO contig00001 1163 523 89 U
    # AAA*GAGAGAAGG*AGGgAGGGaGGAAGgAAGGAAGGAGAG*GAAaGAAA
    # *GAAA*GAGAGAGaGAAAGAAA*GAAA*GG*AAA*GAAGG*AAA*GAAGG
    # ...
    # AF D8QSF4J01CZSRB C 209
    # AF D8QSF4J01CU6TZ U 205
    # AF D8QSF4J01E6204 U 163
    # AF D8QSF4J01A09VN U 214
    # AF D8QSF4J01D2IAU C 176
    # ...
    elsif (/^CO /) {
	# Dump out previous contig
	dump_contig();

	print;

	$cnum++;

	# Skip consensus
	while (<>) {
	    print;
	    last if (/^$/);
	}

	# Skip base quals
	$ _= <>;
	print;
	while (<>) {
	    print;
	    last if (/^$/);
	}

	# Read AF
	while (<>) {
	    last if (/^$/);
	    if (/^AF .*\s(-?[0-9+])/) {
		push(@AF, $_);
	    }
	}
    } elsif (/^BS /) {
	push(@BS, $_);
    } else {
	$seq{$curr_seq} .= $_;
    }
}

dump_contig();

sub dump_contig {
    # Sort AF lines by position
    @AF = map { $_->[0] } sort { $a->[1] <=> $b->[1] }
          map { [$_, / (-?\d+)$/, $_] } @AF;
    
    # Output AF and BS lines
    foreach my $af (@AF) {
	print $af;
    }
    print "\n";
    foreach (@BS) {
	print;
    }

    # Output sequence in order
    print "\n";
    foreach (@AF) {
	/^AF (\S+)/;
	print $seq{$1};
    }

    undef %seq;
    undef @AF;
    undef @BS;
}
