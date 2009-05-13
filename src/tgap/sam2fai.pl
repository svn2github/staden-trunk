#!/usr/bin/perl -w

# Reads a sam file and generates an index suitable for "samtools import".

use strict;

#--- Build index
my %ref; #key=name, value=length
while (<>) {
    my @F=split(" ", $_);
    my $rname = $F[2];

    next if $rname eq "*";
    my $end = $F[3] + length($F[9]);

    $ref{$rname} = $end if (!exists($ref{$rname}) || $ref{$rname} < $end);
}

#--- Output input to stdout
foreach (sort keys %ref) {
    print "$_\t$ref{$_}\n";
}
