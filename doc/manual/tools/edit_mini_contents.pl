#!/usr/bin/perl -w

#
# 14/03/00 jkb
#
# Adds the Home, Up, etc to the mini manual contents pages.
#
while (<>) {
  if (/^<H1>/) {
    print <<EOH
<a href="../index.html"><img src="i/nav_home.gif" alt="home"></a>
<a href="master_contents.html"><img src="i/nav_full.gif" alt="full"></a>
<hr size=4>
EOH
  } elsif (/<H2>Last update on/) {
    next;
  } elsif (/^<HR>$/) {
    print <<EOF
<hr size=4> 
<a href="../index.html"><img src="i/nav_home.gif" alt="home"></a> 
<a href="master_contents.html"><img src="i/nav_full.gif" alt="full"></a> 
EOF
  }
  print;
}
