#!/usr/bin/perl -w

#
# 13/10/95 jkb
#
# Process a list of html files updating any unresolved cross references.
# These may not have been made by texi2html as the reference was in a
# different file.
#
# This assumes the reference will be of the form (on a line by itself)
# "See section `node` in <CITE>section name</CITE>"
# See the maros.texi (fxref, fref) for more details. An assumption has
# also been made that we're using a modified texi2html that writes
# <!-- NODE:name --> and <!-- XREF:node --> lines (which ours does, due to
# mods in texi2html and our macros.texi).
#

%node_table = ();
@FILES = @ARGV;

#
# Loop around all files generating a map of node tables to URLs.
#
while (<ARGV>) {
    if (/^<!-- NODE:/) {
	s/^<!-- NODE://;
	s/ -->\n//;
	$node_name = $_;
    }

    if (/^<H[1-6]><A NAME=/) {
	/NAME="([^"]*)/;
	if ($node_name) {
	    $node_table{$node_name} = "$ARGV#$1";
	    $node_name="";
	}
    }
}

#
# Loop around the files once more updating the cross references
#
while (<@FILES>) {
    print "Processing $_\n";
    $first_line=0;

    $fname = $_;
    open(IN, "< $fname");
    open(OUT, "> _$fname") || die "Couldn't create _$fname";
    while (<IN>) {
        if (/^<!-- XREF:/) {
	    /XREF:(.*) -->$/;
	    $node_name=$1;
	    next;
	}
	if (/section `[^']*' in <CITE>/) {
	    /([^`]*)`([^']*).*<CITE>(.*)<\/CITE>(.*)/;
	    if (exists $node_table{$node_name}) {
	        print "Resolved cross-reference \"$2\"\n";
	        print OUT "$1<A HREF=\"$node_table{$node_name}\">$2</A>$4";
#			  . " in " .
#			  "<A HREF=\"$3_toc.html\">$2_toc.html</A>.";
	        next;
	    } else {
                print "Couldn't resolve \"$2\"\n";
	    }
	}

	s/\n$//;
	if ($first_line) {
	    print OUT "$_";
       	} else {
	    print OUT "\n$_";
	}
    }
    print OUT "\n";
    close(IN);
    close(OUT);

    rename("_$fname", $fname);
}
