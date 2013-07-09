#!/usr/bin/perl

if ($#ARGV < 3) { print "param mismatch\n"; exit 0; }
my $infile=shift;
my $outfile=shift;
my $param=shift;
my $newval=shift;

open(IN,"<$infile") || die("cannot open old runfile for updates\n");
open(OU,">$outfile") || die("cannot open new runfile for updates\n");
while ($line = <IN>) {
    if ($line =~ /^$param=(.*)$/) { print OU "$param=$newval\n"; }
    else { print OU $line; }
}
