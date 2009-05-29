#!/usr/bin/perl -w

# Reads a sam file and generates an index suitable for "samtools import".

use strict;

#--- Build index
my %ref; #key=name, value=length
my %ord;
my $count = 0;
my $last_end = 0;
while (<>) {
    my @F=split(" ", $_);
    my $rname = $F[2];

    next if $rname eq "*";
    my $end = $F[3] + length($F[9]);

    if (!exists($ref{$rname})) {
	$ref{$rname} = $end;
	$ord{$rname} = $count++;
	$last_end = $end;
    } else {
	if ($last_end < $end) {
	    $ref{$rname} = $last_end = $end;
	}
    }
}

#--- Output input to stdout
foreach (sort {$ord{$a} <=> $ord{$b}} keys %ref) {
    print "$_\t$ref{$_}\n";
}
