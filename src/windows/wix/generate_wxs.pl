#!/usr/local/bin/perl

use strict;
use Cwd;
use Win32;

my %shortcuts = ('trev.exe' => 'Trev');

# Component IDs that we need to turn into features.
my @components;

{
    my $globalid=0;
    sub nextid {
        my ($name) = @_;
        $globalid++;
        return "ID$name-$globalid";
   }
}

{
    my $subtime = 0;
    sub guidgen {
        my $time_low = time();
        my $time_mid = $subtime & 0xff;
        my $time_hi_and_version = ($subtime >> 16);
        $time_hi_and_version = ($time_hi_and_version & 0xfff) | 0x4000;
        my $clock_seq_hi_and_reserved = int(rand(0x10000));
        $clock_seq_hi_and_reserved = ($clock_seq_hi_and_reserved & 0x3fff) | 0x8000;
        my $clock_seq_low = int(rand(0x10000));

	$subtime++;

	return sprintf("%08X-%04X-%04X-%04X-0008028F89FB",
	    $time_low, $time_mid, $time_hi_and_version, $clock_seq_hi_and_reserved,
	    $clock_seq_low);
    }
}

sub printheader {
    print "<?xml version=\"1.0\"?>\n";
    print "<Include>\n";
}

sub printdir {
    my ($indent) = @_;
    my $firstfile=1;
    my $sp = "  " x $indent;
    my $cwd = getcwd();

    opendir(my $fh, ".");
    # Add files (in <Component><File/></Component> block
    my @files = readdir($fh);
    foreach (@files) {
	next if /^\.(\.)?$/;
	if (-f "$_") {
	    if ($firstfile) {
		my $id = nextid("dir");
		print "$sp<Component Guid=\"", guidgen(), "\" Id=\"$id\">\n";
		$firstfile = 0;
		push(@components, $id);
	    }
	    my $id = nextid($_);
  	    my $long = "$_";
	    my $short = Win32::GetShortPathName($long);
	    if (exists($shortcuts{$long})) {
	        print "$sp  <File Id=\"$id\" Name=\"$short\" LongName=\"$long\" DiskId=\"1\" src=\"$cwd/$long\">\n";
		print "$sp    <Shortcut Id=\"Shortcut-$long\" Directory=\"StadenMenuFolder\" Name=\"$shortcuts{$long}\" Target=\"[#$id]\" Show=\"minimized\"/>\n";
	        print "$sp  </File>\n";
	    } else {
	        print "$sp  <File Id=\"$id\" Name=\"$short\" LongName=\"$long\" DiskId=\"1\" src=\"$cwd/$long\"/>\n";
	    }
	}
    }
    if (!$firstfile) {
 	print "$sp  <CreateFolder/>\n";
	print "$sp</Component>\n";
    }

    # Recurse down directories.
    foreach (@files) {
	next if /^\.(\.)?$/;
	if (-d "$_") {
	    my $id = nextid($_);
  	    my $long = "$_";
	    my $short = Win32::GetShortPathName($long);
	    print "$sp<Directory Id=\"$id\" Name=\"$short\" LongName=\"$long\">\n";
	    my $lastdir = getcwd();
	    chdir($_);
	    printdir($indent+1);
	    chdir($lastdir);
	    print "$sp</Directory>\n";
	}
    }
    closedir($fh);
}

sub printcomponents {
    print "    <Feature Id=\"main\" Title=\"Staden Package\" Level=\"1\">\n";
    foreach (@components) {
	print "      <ComponentRef Id=\"$_\"/>\n";
    }
    print "    </Feature>\n";
}

sub printfooter {
    print "</Include>\n";
}

printheader();
chdir($ARGV[0]);
print "    <Directory Id=\"TARGETDIR\" Name=\"SourceDir\">\n";
print "      <Directory Id=\"ProgramMenuFolder\" Name=\"PMFolder\">\n";
print "        <Directory Id=\"StadenMenuFolder\" Name=\"Trev\"/>\n";
print "      </Directory>\n";
printdir(2);
print "    </Directory>\n";
printcomponents();
printfooter();
